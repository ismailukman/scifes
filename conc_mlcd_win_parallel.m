clear; clc;

gamma = 1.2;
omega = 1.5;

% ---------- Use all local cores ----------
if isempty(gcp('nocreate'))
    parpool('local');   % start a pool with default (all) workers
end

% ======== INPUTS (two datasets) ========
%mat_path_all     = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_all_groups_1tr_windows.mat';
%mat_path_fes     = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_FES_1tr_windows.mat';
%mat_path_passive = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_Passive_1tr_windows.mat';

mat_path_prefes     = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_fes_1tr_windows.mat';
mat_path_prenfes = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_nfes_1tr_windows.mat';
mat_path_postfes     = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_fes_1tr_windows.mat';
mat_path_postnfes = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_nfes_1tr_windows.mat';


% ======== OUTPUT DIRECTORY ========
out_dir = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/modularity_var/';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% ---------------------------------------------------------
% Helper: run your current analysis pipeline on ONE .mat file
%   (unchanged logic, adds 'tag' to suffix saved variable names)
% ---------------------------------------------------------
function run_one_pipeline(mat_path, out_fullpath, gamma, omega, tag)
    S = load(mat_path);

    % === Gather group variables corr_<group> ===
    fn = fieldnames(S);
    corr_fields = fn(startsWith(fn,'corr_') & ~strcmp(fn,'__meta__'));
    if isempty(corr_fields)
        error('No variables named corr_<group> found in %s', mat_path);
    end

    % Reference region size from the first group
    X0 = S.(corr_fields{1});   % expected (subjects*windows, R, R)
    if ndims(X0) ~= 3
        error('Variable %s must be 3D (subjects*windows, R, R).', corr_fields{1});
    end
    R_ref = size(X0,2);
    if size(X0,3) ~= R_ref
        error('Variable %s has non-square matrices (%dx%d).', corr_fields{1}, R_ref, size(X0,3));
    end

    % Validate and collect in parallel
    stack_list = cell(numel(corr_fields),1);
    parfor i = 1:numel(corr_fields)
        Xi = S.(corr_fields{i});
        if ndims(Xi) ~= 3
            error('Variable %s must be 3D (subjects*windows, R, R).', corr_fields{i});
        end
        if size(Xi,2) ~= R_ref || size(Xi,3) ~= R_ref
            error('Region size mismatch in %s: %dx%d, expected %dx%d.', ...
                  corr_fields{i}, size(Xi,2), size(Xi,3), R_ref, R_ref);
        end
        % Cast to single to reduce memory footprint
        stack_list{i} = single(Xi);    % (Ti, R, R)
    end
    clear X0 Xi

    % Stack along first dim -> (T, R, R), then permute to (R, R, T)
    corr_stack = cat(1, stack_list{:});      % (T, R, R) where T = sum over groups
    clear stack_list S
    corr_all_g = permute(corr_stack, [2 3 1]);  % (R, R, T)
    clear corr_stack

    % ---- Cleaning (vectorized) ----
    % Symmetrize all layers at once
    corr_all_g = 0.5 * (corr_all_g + permute(corr_all_g, [2 1 3]));
    % Set diagonal to 1 for all layers
    N   = size(corr_all_g,1);
    T_g = size(corr_all_g,3);
    maskI = repmat(eye(N,'logical'), 1, 1, T_g);
    corr_all_g(maskI) = 1;
    % Replace non-finite with 0
    corr_all_g(~isfinite(corr_all_g)) = 0;

    % === Convert to multilayer cell array A_g (each cell is R x R) ===
    A_g = squeeze(num2cell(corr_all_g, [1 2]));   %  T_g x 1 cell, each NxN
    % clear corr_all_g   % uncomment to free memory if needed

    % === Multilayer community detection (GENLouvain) ===
    [B_g, twom_g] = multicat(A_g, gamma, omega);
    postprocess_fn = @(Sg) postprocess_categorical_multilayer(Sg, T_g);
    [S_g, Q_g] = iterated_genlouvain(B_g, 10000, 0, 1, 'moverandw', [], postprocess_fn);
    Q_g = Q_g / twom_g;
    S_g = reshape(S_g, N, T_g);
    comm_num_all_g = max(S_g,[],'all');  

    % === Individual multilayer communities ===
    multi_comm_indivi_g = multilayer_community_detection_individual(A_g, 'cat', 100, 'max', gamma, omega);

    % Collect consensus partitions (parallel aggregation)
    n_subj_g = numel(multi_comm_indivi_g);
    N_all_cell = cell(1, n_subj_g);
    parfor subj_i = 1:n_subj_g
        N_all_cell{subj_i} = multi_comm_indivi_g{subj_i}.multi_module_consensus;
    end
    N_all_g = [N_all_cell{:}];

    % === Save results for THIS dataset with suffixed variable names ===
    out = struct();
    out.(sprintf('N_all_g_%s', tag))        = N_all_g;
    out.(sprintf('S_g_%s', tag))            = S_g;
    out.(sprintf('Q_g_%s', tag))            = Q_g;
    out.(sprintf('comm_num_all_g_%s', tag)) = comm_num_all_g;

    save(out_fullpath, '-struct', 'out', '-v7.3');
end


%run_one_pipeline(mat_path_all,     fullfile(out_dir, 'mlcd_all_1tr_wins.mat'),     gamma, omega, 'all');

%run_one_pipeline(mat_path_fes,     fullfile(out_dir, 'mlcd_fes_1tr_wins.mat'),     gamma, omega, 'fes');
%run_one_pipeline(mat_path_passive, fullfile(out_dir, 'mlcd_passive_1tr_wins.mat'), gamma, omega, 'passive');

run_one_pipeline(mat_path_prefes,     fullfile(out_dir, 'mlcd_prefes_1tr_wins.mat'), gamma, omega, 'prefes');
run_one_pipeline(mat_path_prenfes, fullfile(out_dir, 'mlcd_prenfes_1tr_wins.mat'), gamma, omega, 'prenfes');
run_one_pipeline(mat_path_postfes,     fullfile(out_dir, 'mlcd_postfes_1tr_wins.mat'), gamma, omega, 'postfes');
run_one_pipeline(mat_path_postnfes, fullfile(out_dir, 'mlcd_postnfes_1tr_wins.mat'), gamma, omega, 'postnfes');


disp('Processing complete!');
