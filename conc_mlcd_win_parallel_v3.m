% Clear workspace and command window
clear;
clc;

%% Parameter Initialization
gamma = 1.2;
omega = 1.5;

% ---------- Trying to use all local cores ----------
if isempty(gcp('nocreate'))
    parpool('local');   % start a pool with default (all) workers
end

%% MANUAL CONFIGURATION: Number of subjects per group
% Based on my Python output with flatten subjects x windows -> (S*W, R, R) 
n_subjects_g1 = 7;  % PreFES:  (7, 396, 200, 200) -> (2772, 200, 200)
n_subjects_g2 = 5;  % PreNFES: (5, 396, 200, 200) -> (1980, 200, 200)
n_subjects_g3 = 7;  % PostFES: (7, 396, 200, 200) -> (2772, 200, 200)
n_subjects_g4 = 5;  % PostNFES: (5, 396, 200, 200) -> (1980, 200, 200)

fprintf('===== SUBJECT CONFIGURATION =====\n');
fprintf('PreFES:   %d subjects\n', n_subjects_g1);
fprintf('PreNFES:  %d subjects\n', n_subjects_g2);
fprintf('PostFES:  %d subjects\n', n_subjects_g3);
fprintf('PostNFES: %d subjects\n\n', n_subjects_g4);

%% Load adjacency matrices IN PARALLEL
fprintf('Loading data files...\n');
tic;

file_paths = {
    '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_fes_1tr_windows.mat';
    '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_nfes_1tr_windows.mat';
    '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_fes_1tr_windows.mat';
    '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_nfes_1tr_windows.mat'
};

% Parallel load
S_all = cell(4, 1);
corr_fields_all = cell(4, 1);
parfor i = 1:4
    S_all{i} = load(file_paths{i});
    fn = fieldnames(S_all{i});
    corr_fields_all{i} = fn(startsWith(fn,'corr_') & ~strcmp(fn,'__meta__'));
end

fprintf('Data loaded in %.2f seconds\n', toc);

%% Stack and clean correlation matrices IN PARALLEL
fprintf('Processing correlation matrices...\n');
tic;

corr_g_all = cell(4, 1);
T_g_all = zeros(4, 1);

parfor grp = 1:4
    % Each file has ONE variable corr_<group> with shape (subjects*windows, R, R)
    % Just extract the single variable directly
    corr_field = corr_fields_all{grp}{1};  % Get the single corr_ variable name
    corr_data = S_all{grp}.(corr_field);   % (subjects*windows, R, R)
    
    % Permute to (R, R, subjects*windows)
    corr_g = permute(corr_data, [2 3 1]);
    
    % Clean correlation matrices
    corr_g = 0.5 * (corr_g + permute(corr_g, [2 1 3]));
    N = size(corr_g, 1);
    T_g = size(corr_g, 3);
    maskI = repmat(eye(N,'logical'), 1, 1, T_g);
    corr_g(maskI) = 1;
    corr_g(~isfinite(corr_g)) = 0;
    
    corr_g_all{grp} = corr_g;
    T_g_all(grp) = T_g;
end

fprintf('Matrices processed in %.2f seconds\n', toc);

% Extract individual groups
corr_g1 = corr_g_all{1}; T_g1 = T_g_all(1);
corr_g2 = corr_g_all{2}; T_g2 = T_g_all(2);
corr_g3 = corr_g_all{3}; T_g3 = T_g_all(3);
corr_g4 = corr_g_all{4}; T_g4 = T_g_all(4);
N = size(corr_g1, 1);

% Verify subject counts match expectations
fprintf('\nVerifying structure:\n');
fprintf('  PreFES:   %d total windows / %d subjects = %d windows/subject\n', ...
        T_g1, n_subjects_g1, T_g1/n_subjects_g1);
fprintf('  PreNFES:  %d total windows / %d subjects = %d windows/subject\n', ...
        T_g2, n_subjects_g2, T_g2/n_subjects_g2);
fprintf('  PostFES:  %d total windows / %d subjects = %d windows/subject\n', ...
        T_g3, n_subjects_g3, T_g3/n_subjects_g3);
fprintf('  PostNFES: %d total windows / %d subjects = %d windows/subject\n', ...
        T_g4, n_subjects_g4, T_g4/n_subjects_g4);

clear S_all corr_g_all

%% Prepare and reshape adjacency matrices
fprintf('\nPreparing adjacency matrices...\n');
tic;

% Create cell arrays
A_g1 = squeeze(num2cell(corr_g1, [1,2]));
A_g2 = squeeze(num2cell(corr_g2, [1 2]));
A_g3 = squeeze(num2cell(corr_g3, [1 2]));
A_g4 = squeeze(num2cell(corr_g4, [1 2]));

fprintf('\n===== RESHAPING FOR MULTILAYER STRUCTURE =====\n');

% Reshape all groups
A_g1_reshaped = reshape_for_multilayer(A_g1, n_subjects_g1, T_g1, 'PreFES');
A_g2_reshaped = reshape_for_multilayer(A_g2, n_subjects_g2, T_g2, 'PreNFES');
A_g3_reshaped = reshape_for_multilayer(A_g3, n_subjects_g3, T_g3, 'PostFES');
A_g4_reshaped = reshape_for_multilayer(A_g4, n_subjects_g4, T_g4, 'PostNFES');

fprintf('Adjacency prep completed in %.2f seconds\n', toc);

%% Multi-layer Modularity Calculation IN PARALLEL
fprintf('\n===== COMPUTING MULTILAYER MODULARITY =====\n');
tic;

% Store reshaped data in cell array for parallel processing
A_reshaped_all = {A_g1_reshaped, A_g2_reshaped, A_g3_reshaped, A_g4_reshaped};
group_names = {'PreFES', 'PreNFES', 'PostFES', 'PostNFES'};

% Initialize output arrays
S_g_all = cell(4, 1);
Q_g_all = zeros(4, 1);
comm_num_all = zeros(4, 1);

parfor grp = 1:4
    fprintf('Processing %s modularity...\n', group_names{grp});
    
    wins_per_subj = size(A_reshaped_all{grp}, 2);
    
    % Compute modularity for first subject (representative)
    [B_g, twom_g] = multicat(A_reshaped_all{grp}(1,:), gamma, omega);
    postprocess_fn = @(S) postprocess_categorical_multilayer(S, wins_per_subj);
    [S_g, Q_g] = iterated_genlouvain(B_g, 10000, 0, 1, 'moverandw', [], postprocess_fn);
    
    Q_g = Q_g / twom_g;
    S_g = reshape(S_g, N, wins_per_subj);
    comm_num = max(S_g, [], 'all');
    
    S_g_all{grp} = S_g;
    Q_g_all(grp) = Q_g;
    comm_num_all(grp) = comm_num;
end

% Extract individual results
S_g1 = S_g_all{1}; Q_g1 = Q_g_all(1); comm_num_g1 = comm_num_all(1);
S_g2 = S_g_all{2}; Q_g2 = Q_g_all(2); comm_num_g2 = comm_num_all(2);
S_g3 = S_g_all{3}; Q_g3 = Q_g_all(3); comm_num_g3 = comm_num_all(3);
S_g4 = S_g_all{4}; Q_g4 = Q_g_all(4); comm_num_g4 = comm_num_all(4);

fprintf('Modularity computed in %.2f seconds\n', toc);

%% Consensus Community Detection IN PARALLEL
fprintf('\n===== COMPUTING INDIVIDUAL MULTILAYER COMMUNITIES =====\n');
tic;

multi_comm_indivi_all = cell(4, 1);

parfor grp = 1:4
    fprintf('Processing %s individual communities...\n', group_names{grp});
    multi_comm_indivi_all{grp} = multilayer_community_detection_individual(...
        A_reshaped_all{grp}, 'cat', 100, 'max', gamma, omega);
end

% Extract individual results
multi_comm_indivi_g1 = multi_comm_indivi_all{1};
multi_comm_indivi_g2 = multi_comm_indivi_all{2};
multi_comm_indivi_g3 = multi_comm_indivi_all{3};
multi_comm_indivi_g4 = multi_comm_indivi_all{4};

fprintf('Individual communities computed in %.2f seconds\n', toc);

%% Extract consensus partitions IN PARALLEL
fprintf('\nExtracting consensus partitions...\n');
tic;

N_all_g_all = cell(4, 1);

parfor grp = 1:4
    n_subj = numel(multi_comm_indivi_all{grp});
    N_all_cell = cell(1, n_subj);
    
    for subj_i = 1:n_subj
        N_all_cell{subj_i} = multi_comm_indivi_all{grp}{subj_i}.multi_module_consensus;
    end
    
    N_all_g_all{grp} = [N_all_cell{:}];
end

% Extract individual results
N_all_g1 = N_all_g_all{1};
N_all_g2 = N_all_g_all{2};
N_all_g3 = N_all_g_all{3};
N_all_g4 = N_all_g_all{4};

fprintf('Consensus extracted in %.2f seconds\n', toc);

%% Save Results (Sequential)
fprintf('\n===== SAVING RESULTS =====\n');
tic;

out_dir = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/modularity_var/';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Save Group 1 (PreFES)
out_prefes = struct();
out_prefes.N_all_g_prefes = N_all_g1;
out_prefes.S_g_prefes = S_g1;
out_prefes.Q_g_prefes = Q_g1;
out_prefes.comm_num_all_g_prefes = comm_num_g1;
save(fullfile(out_dir, 'mlcd_prefes_1tr_wins.mat'), '-struct', 'out_prefes', '-v7.3');
fprintf('Saved mlcd_prefes_1tr_wins.mat\n');

% Save Group 2 (PreNFES)
out_prenfes = struct();
out_prenfes.N_all_g_prenfes = N_all_g2;
out_prenfes.S_g_prenfes = S_g2;
out_prenfes.Q_g_prenfes = Q_g2;
out_prenfes.comm_num_all_g_prenfes = comm_num_g2;
save(fullfile(out_dir, 'mlcd_prenfes_1tr_wins.mat'), '-struct', 'out_prenfes', '-v7.3');
fprintf('Saved mlcd_prenfes_1tr_wins.mat\n');

% Save Group 3 (PostFES)
out_postfes = struct();
out_postfes.N_all_g_postfes = N_all_g3;
out_postfes.S_g_postfes = S_g3;
out_postfes.Q_g_postfes = Q_g3;
out_postfes.comm_num_all_g_postfes = comm_num_g3;
save(fullfile(out_dir, 'mlcd_postfes_1tr_wins.mat'), '-struct', 'out_postfes', '-v7.3');
fprintf('Saved mlcd_postfes_1tr_wins.mat\n');

% Save Group 4 (PostNFES)
out_postnfes = struct();
out_postnfes.N_all_g_postnfes = N_all_g4;
out_postnfes.S_g_postnfes = S_g4;
out_postnfes.Q_g_postnfes = Q_g4;
out_postnfes.comm_num_all_g_postnfes = comm_num_g4;
save(fullfile(out_dir, 'mlcd_postnfes_1tr_wins.mat'), '-struct', 'out_postnfes', '-v7.3');
fprintf('Saved mlcd_postnfes_1tr_wins.mat\n');

fprintf('Results saved in %.2f seconds\n', toc);


fprintf('\n===== PROCESSING COMPLETE! =====\n');
