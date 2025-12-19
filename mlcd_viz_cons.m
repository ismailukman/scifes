% Load or access the workspace variables
% Assuming multi_comm_indivi_g1, g2, g3, g4 are already loaded in workspace
n_iterations = 100;
n_regions = 200;
n_windows = 50;
group_names = {'cFES-Pre', 'Passive-Pre', 'cFES-Post', 'Passive-Post'};

% Store all group data in a cell array for easy indexing
group_data = {multi_comm_indivi_g1, multi_comm_indivi_g2, multi_comm_indivi_g3, multi_comm_indivi_g4};

% === LOOP THROUGH ALL 4 GROUPS ===
for selected_grp = 1:4
    
    fprintf('\n=== Processing Group: %s ===\n', group_names{selected_grp});
    
    % Preallocate 3D array to store all iterations
    all_modules = zeros(n_regions, n_windows, n_iterations);
    
    % Select the correct group data
    current_group_data = group_data{selected_grp};
    
    % Extract data from each cell (now correctly uses different g1/g2/g3/g4)
    for iter = 1:n_iterations
        all_modules(:, :, iter) = current_group_data{1,1}.multi_module{1, iter};
    end
    
    fprintf('Extracted data shape: %d regions × %d windows × %d iterations\n', ...
            n_regions, n_windows, n_iterations);
    
    % Select 50 evenly-spaced iterations from 100
    n_selected = 50;
    selected_iterations = round(linspace(1, n_iterations, n_selected));
    fprintf('Selected %d iterations evenly from %d total: [%s]\n', ...
            n_selected, n_iterations, num2str(selected_iterations));
    
    % Create figure for 5×10 grid (50 panels for 50 selected iterations)
    fig = figure('Position', [100, 100, 2000, 1000]);
    n_rows = 5;
    n_cols = 10;
    
    % Global colormap scale
    global_min = min(all_modules(:));
    global_max = max(all_modules(:));
    fprintf('Global data range: min = %.2f, max = %.2f\n', global_min, global_max);
    
    % Plot each selected iteration
    for idx = 1:n_selected
        subplot(n_rows, n_cols, idx);
        iter = selected_iterations(idx);
        
        iteration_data = all_modules(:, :, iter);
        
        imagesc(iteration_data);
        clim([global_min, global_max]);
        colormap(parula);
        title(sprintf('Iter %d', iter), 'FontSize', 9, 'FontWeight', 'bold');
        xlabel('Window', 'FontSize', 7);
        ylabel('Region', 'FontSize', 7);
        set(gca, 'FontSize', 7);
        set(gca, 'XTick', [1 10 20 30 40 50]);
        set(gca, 'YTick', [1 50 100 150 200]);
    end
    
    % Single colorbar
    cb = colorbar('Position', [0.92, 0.15, 0.02, 0.7]);
    cb.Label.String = 'Community Assignment';
    cb.Label.FontSize = 12;
    cb.Label.FontWeight = 'bold';
    
    % Main title
    sgtitle(sprintf('%s Subj #1 Community Assignments: 50 Evenly-Selected Iterations - 200 Regions × 50 Windows', ...
            group_names{selected_grp}), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Adjust layout
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.05, 0.05, 0.9, 0.85]);
    set(gcf, 'PaperPositionMode', 'auto');
    
    path_out = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/';
    output_file = fullfile(path_out, sprintf('consensus_iteration_of_community_assignments_viz_%s.png', group_names{selected_grp}));
    
    print(fig, output_file, '-dpng', '-r300');
    fprintf('✓ Saved figure: %s\n', output_file);
    
    set(fig, 'Visible', 'on');
    
end % END GROUP LOOP

% =========================================================================
% Helper function for viridis colormap
% =========================================================================
function cmap = viridis(m)
    if nargin < 1
        m = 256;
    end
    
    viridis_data = [
        0.267004, 0.004874, 0.329415;
        0.282623, 0.140926, 0.457517;
        0.253935, 0.265254, 0.529983;
        0.206756, 0.371758, 0.553117;
        0.163625, 0.471133, 0.558148;
        0.127568, 0.566949, 0.550556;
        0.134692, 0.658636, 0.517649;
        0.266941, 0.748751, 0.440573;
        0.477504, 0.821444, 0.318195;
        0.741388, 0.873449, 0.149561;
        0.993248, 0.906157, 0.143936;
    ];
    
    x = linspace(0, 1, size(viridis_data, 1));
    xi = linspace(0, 1, m);
    cmap = zeros(m, 3);
    
    for i = 1:3
        cmap(:, i) = interp1(x, viridis_data(:, i), xi, 'pchip');
    end
end