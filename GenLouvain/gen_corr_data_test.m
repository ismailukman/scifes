%% Generate REALISTIC Random Correlation Matrices
% Creates brain-like correlation structure to avoid empty communities

clear;
clc;

fprintf('===== GENERATING REALISTIC RANDOM CORRELATION DATA =====\n');

%% Configuration
n_regions = 200;
n_windows = 50;
n_subjects = 1;
n_groups = 4;

group_names = {'pre_fes', 'pre_nfes', 'post_fes', 'post_nfes'};

% Output directory
out_dir = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% Generate Random Data for Each Group
for g = 1:n_groups
    fprintf('\nGenerating data for group: %s\n', group_names{g});
    
    % Initialize correlation matrix: (subjects*windows, R, R)
    corr_data = zeros(n_windows, n_regions, n_regions);
    
    % Define realistic community structure (5 communities)
    n_communities = 5;
    comm_size = n_regions / n_communities;  % 40 regions per community
    community_labels = repelem(1:n_communities, comm_size);
    
    % Generate random correlation matrices for each window
    for w = 1:n_windows
        
        % METHOD: Generate realistic correlation matrix with community structure
        % Step 1: Start with random Gaussian matrix
        A = randn(n_regions, n_regions);
        
        % Step 2: Add community structure (stronger within-community connections)
        for c = 1:n_communities
            idx = find(community_labels == c);
            % Boost within-community connections
            A(idx, idx) = A(idx, idx) + 2.0;  % Stronger signal
        end
        
        % Step 3: Make symmetric
        A = (A + A') / 2;
        
        % Step 4: Compute correlation matrix via covariance-like structure
        % Add small noise to diagonal for numerical stability
        A = A + 0.1 * eye(n_regions);
        
        % Step 5: Convert to correlation (normalize by outer product of std)
        D = diag(1 ./ sqrt(diag(A * A')));
        corr_mat = D * A * D;
        
        % Step 6: Ensure valid correlation properties
        corr_mat = (corr_mat + corr_mat') / 2;  % Perfect symmetry
        corr_mat(1:n_regions+1:end) = 1;        % Diagonal = 1
        
        % Step 7: Clip to valid correlation range [-1, 1]
        corr_mat = max(-0.99, min(0.99, corr_mat));
        
        % Step 8: Add small window-specific variation
        noise = 0.05 * randn(n_regions, n_regions);
        noise = (noise + noise') / 2;
        corr_mat = corr_mat + noise;
        corr_mat = max(-0.99, min(0.99, corr_mat));
        corr_mat(1:n_regions+1:end) = 1;
        
        % Store
        corr_data(w, :, :) = corr_mat;
    end
    
    % Verify shape and properties
    fprintf('  Generated shape: [%d, %d, %d]\n', size(corr_data,1), size(corr_data,2), size(corr_data,3));
    sample_mat = squeeze(corr_data(1,:,:));
    fprintf('  Sample matrix stats: min=%.3f, max=%.3f, mean=%.3f\n', ...
            min(sample_mat(:)), max(sample_mat(:)), mean(sample_mat(:)));
    fprintf('  Diagonal check: all ones? %d\n', all(diag(sample_mat) == 1));
    fprintf('  Symmetry check: max asymmetry = %.6f\n', max(abs(sample_mat - sample_mat'), [], 'all'));
    
    % Create variable name
    safe_g = group_names{g};
    var_name = sprintf('corr_%s', safe_g);
    
    % Save using HDF5 format
    out_file = fullfile(out_dir, sprintf('corr_%s_8tr_windows_1sub_RANDOM.mat', safe_g));
    
    data_to_save = struct();
    data_to_save.(var_name) = corr_data;
    
    save(out_file, '-struct', 'data_to_save', '-v7.3');
    fprintf('  ✓ Saved: %s\n', out_file);
end

fprintf('\n✓ All realistic random data generated successfully!\n');

%% Optional: Visualize one sample to verify structure
fprintf('\n===== VISUALIZATION CHECK =====\n');
test_file = fullfile(out_dir, 'corr_pre_fes_8tr_windows_1sub_RANDOM.mat');
test_data = load(test_file);
field_name = fieldnames(test_data);
test_matrix = test_data.(field_name{1});

% Plot first window
sample_window = squeeze(test_matrix(1, :, :));

figure('Position', [100, 100, 800, 700]);
imagesc(sample_window);
colorbar;
colormap(jet);
clim([-1, 1]);
title('Sample Random Correlation Matrix (Window 1)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Region'); ylabel('Region');
set(gca, 'FontSize', 10);

% Add community boundaries
hold on;
for c = 1:(n_communities-1)
    boundary = c * comm_size + 0.5;
    plot([0.5, n_regions+0.5], [boundary, boundary], 'w--', 'LineWidth', 2);
    plot([boundary, boundary], [0.5, n_regions+0.5], 'w--', 'LineWidth', 2);
end
hold off;

saveas(gcf, fullfile(out_dir, 'random_data_sample.png'));
fprintf('✓ Sample visualization saved: random_data_sample.png\n');

%% Validate that data can be used for community detection
fprintf('\n===== COMMUNITY DETECTION TEST =====\n');
fprintf('Testing if generated data works with genlouvain...\n');

% Quick test on one window
test_window = squeeze(test_matrix(1, :, :));
test_window = (test_window + test_window') / 2;
test_window(~isfinite(test_window)) = 0;

% Compute modularity matrix
k = sum(test_window, 2);
m = sum(k) / 2;
B_test = test_window - (k * k') / (2 * m);

% Run genlouvain
try
    [S_test, Q_test] = genlouvain(B_test, 10000, 0);
    n_comm_test = max(S_test);
    fprintf('✓ Community detection successful!\n');
    fprintf('  Found %d communities, Q = %.4f\n', n_comm_test, Q_test/(2*m));
    fprintf('  Community sizes: %s\n', num2str(histcounts(S_test, 1:n_comm_test+1)));
catch ME
    fprintf('✗ Community detection failed: %s\n', ME.message);
end

fprintf('\n===== DONE =====\n');
fprintf('Random data files ready for use:\n');
for g = 1:n_groups
    fprintf('  %s\n', fullfile(out_dir, sprintf('corr_%s_8tr_windows_1sub_RANDOM.mat', group_names{g})));
end