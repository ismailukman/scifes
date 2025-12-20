% Clear workspace and command window
clear; clc; clear global;


%% Parameter Initialization
gamma = 1.0; %1.2; 1.21; 1
omega = 1.0; %1.5; 0.1; 1

out_dir = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/modularity_var/';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

DO_PLOT = true;   % <-- set to false to skip plotting
density = 0.30;   % keep top 10% of unique (undirected) edges 

% ---------- Trying to use all local cores ----------
if isempty(gcp('nocreate'))
    parpool('local');   % start a pool with default (all) workers
end

%% MANUAL CONFIGURATION: Number of subjects per group
% Based on my Python output with flatten subjects x windows -> (S*W, R, R) 
n_subjects_g1 = 7;  % PreFES:  (7, 50, 200, 200) -> (350, 200, 200)
n_subjects_g2 = 5;  % PreNFES: (5, 50, 200, 200) -> (250, 200, 200)
n_subjects_g3 = 7;  % PostFES: (7, 50, 200, 200) -> (350, 200, 200)
n_subjects_g4 = 5;  % PostNFES: (5, 50, 200, 200) -> (250, 200, 200)

fprintf('===== SUBJECT CONFIGURATION =====\n');
fprintf('PreFES:   %d subjects\n', n_subjects_g1);
fprintf('PreNFES:  %d subjects\n', n_subjects_g2);
fprintf('PostFES:  %d subjects\n', n_subjects_g3);
fprintf('PostNFES: %d subjects\n\n', n_subjects_g4);

%% Load adjacency matrices IN PARALLEL
fprintf('Loading data files...\n');
tic;

file_paths = {
    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_fes_8tr_windows_1sub.mat';
    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_nfes_8tr_windows_1sub.mat';
    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_fes_8tr_windows_1sub.mat';
    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_nfes_8tr_windows_1sub.mat'
    % % 
    '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_fes_8tr_windows.mat';
    '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_nfes_8tr_windows.mat';
    '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_fes_8tr_windows.mat';
    '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_nfes_8tr_windows.mat'

    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_fes_8tr_windows_1sub_RANDOM.mat';
    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_nfes_8tr_windows_1sub_RANDOM.mat';
    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_fes_8tr_windows_1sub_RANDOM.mat';
    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_nfes_8tr_windows_1sub_RANDOM.mat'

    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_fes_8tr_windows_v1.mat';
    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_nfes_8tr_windows_v1.mat';
    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_fes_8tr_windows_v1.mat';
    % '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_nfes_8tr_windows_v1.mat'

};

% Parallel load
S_all = cell(4, 1);
corr_fields_all = cell(4, 1);
parfor i = 1:4
    S_all{i} = load(file_paths{i});
    fn = fieldnames(S_all{i});
    corr_fields_all{i} = fn(startsWith(fn,'corr_') & ~strcmp(fn,'__meta__'));
end

%% Stack and clean correlation matrices IN PARALLEL
fprintf('Processing correlation matrices...\n');

corr_g_all = cell(4, 1);
T_g_all = zeros(4, 1);
adj_g_all = cell(4, 1);
n_subjects_array = [n_subjects_g1, n_subjects_g2, n_subjects_g3, n_subjects_g4];

parfor grp = 1:4
% for grp = 1:4
    % Just extract the single variable directly
    corr_field = corr_fields_all{grp}{1};      % single corr_ variable name
    corr_data  = S_all{grp}.(corr_field);      % (subjects*windows, R, R)

    % Permute to (R, R, subjects*windows)
    corr_g = permute(corr_data, [2 3 1]);     

    % Step 1: Symmetry correction
    asym_before = corr_g - permute(corr_g, [2 1 3]);
    max_asym    = max(abs(asym_before), [], 'all');
    n_asym      = nnz(abs(asym_before) > eps);
    corr_g      = 0.5 * (corr_g + permute(corr_g,[2 1 3]));

    % Step 2: Diagonal replacement
    N   = size(corr_g,1);
    T_g = size(corr_g,3);
    maskI = repmat(eye(N,'logical'),1,1,T_g);
    n_diag_off = nnz(corr_g(maskI) ~= 1);
    corr_g(maskI) = 0;   % (your choice: 0 here since you remove self-edges later)

    % Step 3: Non-finite cleanup
    n_nonfinite = nnz(~isfinite(corr_g));
    corr_g(~isfinite(corr_g)) = 0;

    % Report summary
    fprintf('Cleaning summary (group %d):\n', grp);
    fprintf('  Symmetry: max asym = %.3e, mismatched entries = %d\n', max_asym, n_asym);
    fprintf('  Diagonal corrections: %d entries set to 0\n', n_diag_off);
    fprintf('  Non-finite corrections: %d entries set to 0\n', n_nonfinite);

    %% === Build adjacency matrices (proportional thresholding, density; weighted, no binarization) ===
    adj_g   = zeros(N, N, T_g);

    upmask = triu(true(N), 1);               % upper triangle (no diag)
    [ri, ci] = find(upmask);
    M = numel(ri);                           % # undirected edges
    K = max(0, round(density * M));          % target edges per window

    for t = 1:T_g
        Ak = corr_g(:,:,t);
        Ak(1:N+1:end) = 0;                   % no self-edges
        Ak(Ak < 0)    = 0;                   % droe  negatives for unsigned

        % Gather upper-tri weights, keep the top K positives
        w = Ak(upmask);
        pos_idx = find(w > 0);
        if isempty(pos_idx) || K == 0
            adj_g(:,:,t) = zeros(N);
            continue;
        end
        Kt = min(K, numel(pos_idx));
        [~, ord] = sort(w(pos_idx), 'descend');
        keep_pos = pos_idx(ord(1:Kt));

        Akeep = zeros(N);
        lin   = sub2ind([N N], ri(keep_pos), ci(keep_pos));
        Akeep(lin) = w(keep_pos);
        Akeep = Akeep + Akeep.';             % symmetrize
        adj_g(:,:,t) = Akeep;
    end

    % achieved density report (upper-tri)
    upmask3 = repmat(upmask, 1, 1, T_g);
    nz_upper = squeeze(sum(sum((adj_g > 0) & upmask3, 1), 2));  % 1×T_g
    achieved_density = nz_upper / M;
    fprintf('Adjacency build (group %d): target=%.3f, median=%.3f (IQR: [%.3f, %.3f])\n', ...
        grp, density, median(achieved_density), quantile(achieved_density,0.25), quantile(achieved_density,0.75));

    % Combined visualization: Raw vs Cleaned vs Adjacency
    if DO_PLOT
        % Raw stack (R x R x T_total)
        sz = size(corr_data);
        if numel(sz) ~= 3
            error('corr_data must be 3-D (either T×R×R or R×R×T).');
        end
        if sz(1) == sz(2)
            corr_stack = corr_data;          % already R x R x T
            T_total    = sz(3);
        else
            corr_stack = permute(corr_data, [2 3 1]);   % R x R x T
            T_total    = sz(1);
        end

        % Row count depends on whether adjacency is available
        show_adj = ~isempty(adj_g);
        nrows    = 2 + show_adj;             % 2 or 3 rows

        % Indices to plot: 10, 20, 30, 40, 50 (bounded by T)
        T_use = min([T_total, size(corr_g,3), size(adj_g,3)]);
        idx   = 10:10:min(50, T_use);
        if isempty(idx)
            warning('No indices selected (T = %d). Need at least 10 windows.', T_use);
        else
            % Shared color limits across selected panels
            min_all = inf; max_all = -inf;
            for j = 1:numel(idx)
                kWin = idx(j);
                Araw = corr_stack(:,:,kWin);
                Acln = corr_g(:,:,kWin);
                min_all = min([min_all, min(Araw(:)), min(Acln(:))]);
                max_all = max([max_all, max(Araw(:)), max(Acln(:))]);
                if show_adj
                    Aadj = adj_g(:,:,kWin);
                    min_all = min(min_all, min(Aadj(:)));
                    max_all = max(max_all, max(Aadj(:)));
                end
            end
            clims = [min_all, max_all];

            % Layout: rows = Raw; Cleaned; (Adjacency), columns = #windows
            nplots = numel(idx);
            ncols  = nplots;

            % Create invisible figure (then save)
            f  = figure('Visible','off','Color','w', ...
                'Name', sprintf('Group %d: Raw vs Cleaned vs Adjacency', grp));
            tl = tiledlayout(nrows, ncols, 'Padding','compact', 'TileSpacing','compact');

            for j = 1:nplots
                kWin = idx(j);

                % Row 1: RAW
                ax1 = nexttile(j);
                imagesc(corr_stack(:,:,kWin), clims);
                axis(ax1, 'image'); axis(ax1, 'off');
                title(ax1, sprintf('Raw %d (min=%.2f, max=%.2f)', kWin, ...
                    min(corr_stack(:,:,kWin),[],'all'), max(corr_stack(:,:,kWin),[],'all')), 'FontSize', 6);

                % Row 2: CLEANED
                ax2 = nexttile(j + ncols);
                imagesc(corr_g(:,:,kWin), clims);
                axis(ax2, 'image'); axis(ax2, 'off');
                title(ax2, sprintf('Clean %d (min=%.2f, max=%.2f)', kWin, ...
                    min(corr_g(:,:,kWin),[],'all'), max(corr_g(:,:,kWin),[],'all')), 'FontSize', 6);

                % Row 3: ADJACENCY
                if show_adj
                    ax3 = nexttile(j + 2*ncols);
                    imagesc(adj_g(:,:,kWin), clims);
                    axis(ax3, 'image'); axis(ax3, 'off');
                    title(ax3, sprintf('Adj %d (min=%.2f, max=%.2f)', kWin, ...
                        min(adj_g(:,:,kWin),[],'all'), max(adj_g(:,:,kWin),[],'all')), 'FontSize', 6);
                end
            end

            colormap(parula);
            cb = colorbar;            % shared colorbar on the right for tiledlayout 
            try, cb.Layout.Tile = 'east'; end

            if show_adj
                title(tl, 'FC — Raw (top) vs Cleaned (middle) vs Adjacency (bottom);', 'FontSize', 7);
            else
                title(tl, 'FC — Raw (top) vs Cleaned (bottom); shared color scale', 'FontSize', 7);
            end

            % Save directly to file 
            out_plot_dir = fullfile(out_dir, 'plots');
            if ~exist(out_plot_dir, 'dir'), mkdir(out_plot_dir); end
            out_file = fullfile(out_plot_dir, sprintf('group%d_plot.png', grp));
            exportgraphics(f, out_file, 'Resolution', 200);

            close(f);  % free memory
        end
    end

    % Store per-group results  
    corr_g_all{grp} = corr_g;
    T_g_all(grp)    = T_g;
    adj_g_all{grp}  = adj_g;
end




% Extract individual groups
corr_g1 = corr_g_all{1}; T_g1 = T_g_all(1); adj_g1 = adj_g_all{1};
corr_g2 = corr_g_all{2}; T_g2 = T_g_all(2); adj_g2 = adj_g_all{2};
corr_g3 = corr_g_all{3}; T_g3 = T_g_all(3); adj_g3 = adj_g_all{3};
corr_g4 = corr_g_all{4}; T_g4 = T_g_all(4); adj_g4 = adj_g_all{4};
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

% Create cell arrays
A_g1 = squeeze(num2cell(adj_g1, [1,2]));
A_g2 = squeeze(num2cell(adj_g2, [1 2]));
A_g3 = squeeze(num2cell(adj_g3, [1 2]));
A_g4 = squeeze(num2cell(adj_g4, [1 2]));

fprintf('\n===== RESHAPING FOR MULTILAYER STRUCTURE =====\n');

% Reshape all groups
A_g1_reshaped = reshape_for_multilayer(A_g1, n_subjects_g1, T_g1, 'PreFES');
A_g2_reshaped = reshape_for_multilayer(A_g2, n_subjects_g2, T_g2, 'PreNFES');
A_g3_reshaped = reshape_for_multilayer(A_g3, n_subjects_g3, T_g3, 'PostFES');
A_g4_reshaped = reshape_for_multilayer(A_g4, n_subjects_g4, T_g4, 'PostNFES');


%% Multi-layer Modularity Calculation IN PARALLEL
fprintf('\n===== COMPUTING MULTILAYER MODULARITY =====\n');

% Store reshaped data in cell array for parallel processing
A_reshaped_all = {A_g1_reshaped, A_g2_reshaped, A_g3_reshaped, A_g4_reshaped};
n_subjects_all = [n_subjects_g1, n_subjects_g2, n_subjects_g3, n_subjects_g4];
group_names = {'PreFES', 'PreNFES', 'PostFES', 'PostNFES'};

% Initialize output arrays
% S_g_all = cell(4, 1);
% Q_g_all = zeros(4, 1);
comm_num_all = zeros(4, 1);

% parfor grp = 1:4
%     fprintf('Processing %s modularity...\n', group_names{grp});
% 
%     n_subj = n_subjects_all(grp);
%     wins_per_subj = size(A_reshaped_all{grp}, 2);
% 
%     % Initialize S_g to collect all subjects
%     S_g_full = [];
% 
%     % Compute modularity for ALL subjects in this group
%     for subj_i = 1:n_subj
%         % [B_g, twom_g] = multicat(A_reshaped_all{grp}(subj_i,:), gamma, omega);
%         [B_g, twom_g] = multiord(A_reshaped_all{grp}(subj_i,:), gamma, omega);
%         % postprocess_fn = @(S) postprocess_categorical_multilayer(S, wins_per_subj);
%         postprocess_fn = @(S) postprocess_ordinal_multilayer(S, wins_per_subj);
%         [S_g, Q_g] = iterated_genlouvain(B_g, 10000, 0, 1, 'moverandw', [], postprocess_fn);
% 
%         S_g = reshape(S_g, N, wins_per_subj);
%         S_g_full = [S_g_full, S_g];  % Concatenate across subjects
%     end
% 
%     % Use first subject's metrics as representative
%     % [B_g, twom_g] = multicat(A_reshaped_all{grp}(1,:), gamma, omega);
%     [B_g, twom_g] = multiord(A_reshaped_all{grp}(1,:), gamma, omega);
%     % postprocess_fn = @(S) postprocess_categorical_multilayer(S, wins_per_subj);
%     postprocess_fn = @(S) postprocess_ordinal_multilayer(S, wins_per_subj);
%     [~, Q_g] = iterated_genlouvain(B_g, 10000, 0, 1, 'moverandw', [], postprocess_fn);
%     Q_g = Q_g / twom_g;
% 
%     comm_num = max(S_g_full, [], 'all');
% 
%     S_g_all{grp} = S_g_full;
%     Q_g_all(grp) = Q_g;
%     comm_all(grp) = comm_num;
% end
% 
% % Extract individual results
% S_g1 = S_g_all{1}; Q_g1 = Q_g_all(1); comm_g1 = comm_all(1);
% S_g2 = S_g_all{2}; Q_g2 = Q_g_all(2); comm_g2 = comm_all(2);
% S_g3 = S_g_all{3}; Q_g3 = Q_g_all(3); comm_g3 = comm_all(3);
% S_g4 = S_g_all{4}; Q_g4 = Q_g_all(4); comm_g4 = comm_all(4);


%% Consensus Community Detection IN PARALLEL
fprintf('\n===== COMPUTING INDIVIDUAL MULTILAYER COMMUNITIES =====\n');

multi_comm_indivi_all = cell(4, 1);

parfor grp = 1:4
    fprintf('Processing %s individual communities...\n', group_names{grp});
    multi_comm_indivi_all{grp} = multilayer_community_detection_individual( ...
    A_reshaped_all{grp}, 'ord', 'n_repeat', 100, 'thresh_type', 'max', 'gamma', gamma, 'omega', omega); % cat

end

% Extract individual results
multi_comm_indivi_g1 = multi_comm_indivi_all{1};
multi_comm_indivi_g2 = multi_comm_indivi_all{2};
multi_comm_indivi_g3 = multi_comm_indivi_all{3};
multi_comm_indivi_g4 = multi_comm_indivi_all{4};


%% Extract consensus partitions IN PARALLEL
fprintf('\nExtracting consensus partitions...\n');

N_all_g_all = cell(4, 1);
comm_num_g_all = cell(4, 1);
Qmod_g_all = cell(4, 1);

parfor grp = 1:4
    n_subj = numel(multi_comm_indivi_all{grp});
    N_all_cell = cell(1, n_subj);
    comm_num_cell = cell(1, n_subj); Qmod = cell(1, n_subj)
    
    for subj_i = 1:n_subj
        subj_data = multi_comm_indivi_all{grp}{subj_i};
        N_all_cell{subj_i} = subj_data.multi_module_consensus;
        comm_num_cell{subj_i} = max(subj_data.multi_comm_consensus);
        Qmod{subj_i} = mode(cell2mat(subj_data.multi_modQ));
        
    end
    
    N_all_g_all{grp} = [N_all_cell{:}];
    comm_num_g_all{grp} = [comm_num_cell{:}];
    Qmod_g_all{grp} = [Qmod{:}];
end

% Extract individual group results
N_all_g1 = N_all_g_all{1};
N_all_g2 = N_all_g_all{2};
N_all_g3 = N_all_g_all{3};
N_all_g4 = N_all_g_all{4};

comm_cons_g1 = comm_num_g_all{1};
comm_cons_g2 = comm_num_g_all{2};
comm_cons_g3 = comm_num_g_all{3};
comm_cons_g4 = comm_num_g_all{4};

Qmod_g1 = Qmod_g_all{1};
Qmod_g2 = Qmod_g_all{2};
Qmod_g3 = Qmod_g_all{3};
Qmod_g4 = Qmod_g_all{4};

% %% Verify shapes match
% fprintf('\n===== SHAPE VERIFICATION =====\n');
% fprintf('S_g1: %d x %d, N_all_g1: %d x %d\n', size(S_g1,1), size(S_g1,2), size(N_all_g1,1), size(N_all_g1,2));
% fprintf('S_g2: %d x %d, N_all_g2: %d x %d\n', size(S_g2,1), size(S_g2,2), size(N_all_g2,1), size(N_all_g2,2));
% fprintf('S_g3: %d x %d, N_all_g3: %d x %d\n', size(S_g3,1), size(S_g3,2), size(N_all_g3,1), size(N_all_g3,2));
% fprintf('S_g4: %d x %d, N_all_g4: %d x %d\n', size(S_g4,1), size(S_g4,2), size(N_all_g4,1), size(N_all_g4,2));

%% Save Results (Sequential)
fprintf('\n===== SAVING RESULTS =====\n');



% Save Group 1 (PreFES)
out_prefes = struct();
out_prefes.N_all_g_prefes = N_all_g1;
% out_prefes.S_g_prefes = S_g1;
out_prefes.Q_g_prefes = Qmod_g1;
out_prefes.comm_cons_all_g_prefes = comm_cons_g1;
save(fullfile(out_dir, 'mlcd_prefes_8tr_wins.mat'), '-struct', 'out_prefes', '-v7.3');
fprintf('Saved mlcd_prefes_8tr_wins.mat\n');

% Save Group 2 (PreNFES)
out_prenfes = struct();
out_prenfes.N_all_g_prenfes = N_all_g2;
% out_prenfes.S_g_prenfes = S_g2;
out_prenfes.Q_g_prenfes = Qmod_g2;
out_prenfes.comm_cons_all_g_prenfes = comm_cons_g2;
save(fullfile(out_dir, 'mlcd_prenfes_8tr_wins.mat'), '-struct', 'out_prenfes', '-v7.3');
fprintf('Saved mlcd_prenfes_8tr_wins.mat\n');

% Save Group 3 (PostFES)
out_postfes = struct();
out_postfes.N_all_g_postfes = N_all_g3;
% out_postfes.S_g_postfes = S_g3;
out_postfes.Q_g_postfes = Qmod_g3;
out_postfes.comm_cons_all_g_postfes = comm_cons_g3;
save(fullfile(out_dir, 'mlcd_postfes_8tr_wins.mat'), '-struct', 'out_postfes', '-v7.3');
fprintf('Saved mlcd_postfes_8tr_wins.mat\n');

% Save Group 4 (PostNFES)
out_postnfes = struct();
out_postnfes.N_all_g_postnfes = N_all_g4;
% out_postnfes.S_g_postnfes = S_g4;
out_postnfes.Q_g_postnfes = Qmod_g4;
out_postnfes.comm_cons_all_g_postnfes = comm_cons_g4;
save(fullfile(out_dir, 'mlcd_postnfes_8tr_wins.mat'), '-struct', 'out_postnfes', '-v7.3');
fprintf('Saved mlcd_postnfes_8tr_wins.mat\n');

fprintf('Done in %.2f seconds\n', toc);

