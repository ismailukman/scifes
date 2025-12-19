% Load or access the workspace variables
% Assuming S_g1, S_g2 are already loaded in workspace

n_regions = 200;
n_subjects_vec = [50, 50, 50, 50];  % 
% n_subjects_vec = [1, 1, 1, 1];  % 

% Group names corresponding to S_g1, S_g2
group_names = {'cFES-Pre', 'Passive-Pre', 'cFES-Post', 'Passive-Post'};
type_id = {'consensus', 'non_consensus'};
% group_vars = {S_g1, S_g2, S_g3, S_g4};
group_vars = {multi_comm_indivi_g1{1,1}.multi_module_consensus, multi_comm_indivi_g2{1,1}.multi_module_consensus, multi_comm_indivi_g3{1,1}.multi_module_consensus, multi_comm_indivi_g4{1,1}.multi_module_consensus};


% Validate shapes
for g = 1:4
    [r, c] = size(group_vars{g});
    fprintf('S_g%d (%s): %d regions × %d subjects\n', g, group_names{g}, r, c);
    assert(r == n_regions && c == n_subjects_vec(g), ...
           'Expected shape 200×%d, got %d×%d for %s', n_subjects_vec(g), r, c, group_names{g});
end

% Process each group (single plot per figure)
for g = 1:4
    % Extract data for this group (200 regions × 22/23 subjects)
    group_data = group_vars{g};
    group_name = group_names{g};
    n_subjects = size(group_data, 2);  % actual number of subjects for this group
    
    % Create figure; single panel
    fig = figure('Position', [100, 100, 2000, 1000]);
    
    % Track global min/max for consistent colormap
    global_min = min(group_data(:));
    global_max = max(group_data(:));
    
    fprintf('\n%s: Global data range: min = %.2f, max = %.2f\n', ...
            group_name, global_min, global_max);
    
    % Single heatmap of region (rows) × subject (columns)
    imagesc(group_data);
    clim([global_min, global_max]);
    
    % Colormap (choose one)
    % colormap(viridis);
    colormap(parula);
    % colormap(turbo);
    % colormap(jet);
    
    % Labels and title
    % title(sprintf('%s: Community Label Assignments', group_name), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Windows', 'FontSize', 12);
    ylabel('Region',  'FontSize', 12);
    set(gca, 'FontSize', 10);
    
    % Ticks
    if n_subjects >= 6
        set(gca, 'XTick', unique([1, round(n_subjects/4), round(n_subjects/2), round(3*n_subjects/4), n_subjects]));
    else
        set(gca, 'XTick', 1:n_subjects);
    end
    set(gca, 'YTick', [1, 50, 100, 150, 200]);
    
    % Single colorbar
    cb = colorbar('Position', [0.92, 0.15, 0.02, 0.7]);
    cb.Label.String = 'Community Assignment';
    cb.Label.FontSize = 12;
    cb.Label.FontWeight = 'bold';
    
    % Main title with shape info
    sgtitle(sprintf('%s - Community Label Assignments 200 Regions × %d Windows', ...
            group_name, n_subjects), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Layout and save
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.05, 0.05, 0.9, 0.85]);
    set(gcf, 'PaperPositionMode', 'auto');

    path_out = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/';
    output_file = fullfile(path_out, sprintf('community_assignments_%s_viz_%s.png', group_name, type_id{1}));
    
    % output_file = sprintf('community_assignments_%s_viz_%s.png', group_name, type_id{1});
    
    print(fig, output_file, '-dpng', '-r300');
    fprintf('✓ Figure saved as: %s\n', output_file);
    
    set(fig, 'Visible', 'on');
end

fprintf('\n✓ Both groups processed successfully!\n');

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
    x  = linspace(0, 1, size(viridis_data, 1));
    xi = linspace(0, 1, m);
    cmap = zeros(m, 3);
    for i = 1:3
        cmap(:, i) = interp1(x, viridis_data(:, i), xi, 'pchip');
    end
end
