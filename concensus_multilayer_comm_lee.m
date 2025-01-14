% Clear workspace and command window
clear;
clc;

%% Parameter Initialization
gamma = 1.21; % Resolution parameter
log_omega = -1; % Logarithm of inter-layer coupling
omega = 10^log_omega; % Inter-layer coupling value

% Load adjacency matrices 
corr_g1 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_fes.mat'))), [2 3 1]);
corr_g2 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_fes.mat'))), [2 3 1]);

% Group PreFES
% Prepare adjacency/connectivity matrices
A_g1 = squeeze(num2cell(corr_g1, [1, 2])); % Cell array of adjacency matrices
A_g2 = squeeze(num2cell(corr_g2,[1 2]));  

N = length(A_g1{1}); % Number of nodes
T_g1 = size(corr_g1,3); % Number of layers (subjects)
T_g2 = size(corr_g2,3);

%% Multi-layer Modularity Calculation
[B_g1, twom_g1] = multicat(A_g1, gamma, omega);
postprocess_fn1 = @(S_g1) postprocess_categorical_multilayer(S_g1, T_g1);
% Run multi-layer modularity optimization
[S_g1, Q_g1] = iterated_genlouvain(B_g1, 10000, 0, 1, 'moverandw', [], postprocess_fn1);
Q_g1 = Q_g1 / twom_g1; % Normalize modularity score
S_g1 = reshape(S_g1, N, T_g1); % Reshape results
% Extract community structure
comm_num_g1 = max(S_g1, [], 'all'); % Number of communities

% Group PostFES
[B_g2,twom_g2] = multicat(A_g2, gamma, omega);
postprocess_fn2 = @(S_g2)postprocess_categorical_multilayer(S_g2, T_g2);
[S_g2, Q_g2] = iterated_genlouvain(B_g2, 10000, 0, 1, 'moverandw', [], postprocess_fn2);
Q_g2 = Q_g2/twom_g2;
S_g2 = reshape(S_g2, N, T_g2);
comm_num_g2 = max(S_g2,[],'all'); % number of communities

%% Consensus Community Detection
% multi-layer community structure within individuals (individual level)
multi_comm_indivi_g1 = multilayer_community_detection_individual(...
    A_g1, 'cat', 100, 'max', gamma, omega);  % A_g1, 'cat', 100, 'max', gamma, omega

multi_comm_indivi_g2 = multilayer_community_detection_individual(...
    A_g2, 'cat', 100, 'max', gamma, omega); 

N_all_g1 = [];
N_all_g2 = [];
n_subj_g1 = numel(multi_comm_indivi_g1);
n_subj_g2 = numel(multi_comm_indivi_g2);
% n_node = size(multi_comm_indivi_g1{1}.multi_module_consensus, 1);

for subj_i = 1:n_subj_g1
    % Extract and append consensus community structure for each subject
    N_all_g1 = [N_all_g1, multi_comm_indivi_g1{subj_i}.multi_module_consensus]; 
end

for subj_i = 1:n_subj_g2
    % Extract and append consensus community structure for each subject
    N_all_g2 = [N_all_g2, multi_comm_indivi_g2{subj_i}.multi_module_consensus]; 
end

% multi_module_consensus_sorted

% pause(2);
% consensus community structure across individuals (group level)
multi_comm_group_g1 = multilayer_community_detection_group(multi_comm_indivi_g1, ...
    100, 'max', gamma);
multi_comm_group_g2 = multilayer_community_detection_group(multi_comm_indivi_g2, ...
    100, 'max', gamma);
N_all_group_g1 = multi_comm_group_g1.multi_module_consensus;
N_all_group_g2 = multi_comm_group_g2.multi_module_consensus;

%% Save Results
% output_dir = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/modularity_var/';

cd '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/modularity_var/';
 
filename = sprintf('S_indivi_g1_gamma%.1f,%.1f.mat', gamma, log_omega); save(filename,'N_all_g1');
filename = sprintf('S_group_g1_gamma%.1f,%.1f.mat', gamma, log_omega); save(filename,'N_all_group_g1');

filename = sprintf('S_indivi_g2_gamma%.1f,%.1f.mat', gamma, log_omega); save(filename,'N_all_g2');
filename = sprintf('S_group_g2_gamma%.1f,%.1f.mat', gamma, log_omega); save(filename,'N_all_group_g2');

filename = sprintf('S_prefes_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g1'); 
filename = sprintf('S_postfes_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g2');


disp('Processing complete!');
