% Clear workspace and command window
clear;
clc;

%% 7a.Parameter Initialization
% gamma = 1.21; % Resolution parameter
% log_omega = -1; % Logarithm of inter-layer coupling
% omega = 10^log_omega; % Inter-layer coupling value

gamma = 1.2;
omega = 1.5;

%% 7b. Load adjacency matrices 
corr_g1 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_fes.mat'))), [2 3 1]);
corr_g2 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_nfes.mat'))), [2 3 1]);
corr_g3 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_fes.mat'))), [2 3 1]); 
corr_g4 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_nfes.mat'))), [2 3 1]);

% All group
corr_all_g = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_all_groups.mat'))), [2 3 1]);

% Group Pre-Post | FES-c
% Prepare adjacency/connectivity matrices
A_g1 = squeeze(num2cell(corr_g1, [1,2])); % Cell array of adjacency matrices

A_g2 = squeeze(num2cell(corr_g2, [1 2])); 
A_g3 = squeeze(num2cell(corr_g3, [1 2]));  
A_g4 = squeeze(num2cell(corr_g4, [1 2]));  

A_g = squeeze(num2cell(corr_all_g, [1 2])); 

N = length(A_g1{1}); % Number of nodes
T_g1 = size(corr_g1,3); % Number of layers (subjects)
T_g2 = size(corr_g2,3);
T_g3 = size(corr_g3,3);
T_g4 = size(corr_g4,3);
T_g = size(corr_all_g,3);

%% 7c. Multi-layer Modularity Calculation (Our graph theory method) 
% Group PreFES
[B_g1, twom_g1] = multicat(A_g1, gamma, omega);
postprocess_fn1 = @(S_g1) postprocess_categorical_multilayer(S_g1, T_g1);
% Run multi-layer modularity optimization
[S_g1, Q_g1] = iterated_genlouvain(B_g1, 10000, 0, 1, 'moverandw', [], postprocess_fn1);
Q_g1 = Q_g1 / twom_g1; % Normalize modularity score
S_g1 = reshape(S_g1, N, T_g1); % Reshape results
% Extract community structure
comm_num_g1 = max(S_g1, [], 'all'); % Number of communities

% Group PreSham
[B_g2,twom_g2] = multicat(A_g2, gamma, omega);
postprocess_fn2 = @(S_g2)postprocess_categorical_multilayer(S_g2, T_g2);
[S_g2, Q_g2] = iterated_genlouvain(B_g2, 10000, 0, 1, 'moverandw', [], postprocess_fn2);
Q_g2 = Q_g2/twom_g2;
S_g2 = reshape(S_g2, N, T_g2);
comm_num_g2 = max(S_g2,[],'all'); % number of communities

% Group PostFES
[B_g3,twom_g3] = multicat(A_g3, gamma, omega);
postprocess_fn3 = @(S_g3)postprocess_categorical_multilayer(S_g3, T_g3);
[S_g3, Q_g3] = iterated_genlouvain(B_g3, 10000, 0, 1, 'moverandw', [], postprocess_fn3);
Q_g3 = Q_g3/twom_g3;
S_g3 = reshape(S_g3, N, T_g3);
comm_num_g3 = max(S_g3,[],'all'); % number of communities

% Group PostSham
[B_g4,twom_g4] = multicat(A_g4, gamma, omega);
postprocess_fn4 = @(S_g4)postprocess_categorical_multilayer(S_g4, T_g4);
[S_g4, Q_g4] = iterated_genlouvain(B_g4, 10000, 0, 1, 'moverandw', [], postprocess_fn4);
Q_g4 = Q_g4/twom_g4;
S_g4 = reshape(S_g4, N, T_g4);
comm_num_g4 = max(S_g4,[],'all'); % number of communities

% Group All 
[B_g, twom_g] = multicat(A_g, gamma, omega);
postprocess_fn = @(S_g)postprocess_categorical_multilayer(S_g, T_g);
[S_g, Q_g] = iterated_genlouvain(B_g, 10000, 0, 1, 'moverandw', [], postprocess_fn);
Q_g = Q_g/twom_g;
S_g = reshape(S_g, N, T_g);
comm_num_all_g = max(S_g,[],'all'); % number of communities

%% 7d. Consensus Community Detection (Lee et al method)
%% 7di. multi-layer community structure within individuals (individual level)
% multi_comm_indivi_g1 = multilayer_community_detection_individual(...
%     A_g1, 'cat', 100, 'max', gamma, omega);  % A_g1, 'cat', 100, 'max', gamma, omega
% 
% multi_comm_indivi_g2 = multilayer_community_detection_individual(...
%     A_g2, 'cat', 100, 'max', gamma, omega); 
% 
% 
% multi_comm_indivi_g3 = multilayer_community_detection_individual(...
%     A_g3, 'cat', 100, 'max', gamma, omega); 
% 
% multi_comm_indivi_g4 = multilayer_community_detection_individual(...
%     A_g4, 'cat', 100, 'max', gamma, omega); 

multi_comm_indivi_g = multilayer_community_detection_individual(A_g, 'cat', 100, 'max', gamma, omega); 

% N_all_g1 = [];
% N_all_g2 = [];
% N_all_g3 = [];
% N_all_g4 = [];
N_all_g = [];
% n_subj_g1 = numel(multi_comm_indivi_g1);
% n_subj_g2 = numel(multi_comm_indivi_g2);
% n_subj_g3 = numel(multi_comm_indivi_g3);
% n_subj_g4 = numel(multi_comm_indivi_g4);
n_subj_g = numel(multi_comm_indivi_g);
% n_node = size(multi_comm_indivi_g1{1}.multi_module_consensus, 1);

% for subj_i = 1:n_subj_g1
%     % Extract and append consensus community structure for each subject
%     N_all_g1 = [N_all_g1, multi_comm_indivi_g1{subj_i}.multi_module_consensus]; 
% end
% 
% for subj_i = 1:n_subj_g2
%     % Extract and append consensus community structure for each subject
%     N_all_g2 = [N_all_g2, multi_comm_indivi_g2{subj_i}.multi_module_consensus]; 
% end
% 
% for subj_i = 1:n_subj_g3
%     % Extract and append consensus community structure for each subject
%     N_all_g3 = [N_all_g3, multi_comm_indivi_g3{subj_i}.multi_module_consensus]; 
% end
% 
% for subj_i = 1:n_subj_g4
%     % Extract and append consensus community structure for each subject
%     N_all_g4 = [N_all_g4, multi_comm_indivi_g4{subj_i}.multi_module_consensus]; 
% end

for subj_i = 1:n_subj_g
    % Extract and append consensus community structure for each subject
    N_all_g = [N_all_g, multi_comm_indivi_g{subj_i}.multi_module_consensus]; 
end
% multi_module_consensus_sorted

% pause(2);
%% 7dii. Consensus community structure across individuals (group level)
% multi_comm_group_g1 = multilayer_community_detection_group(multi_comm_indivi_g1, ...
%     100, 'max', gamma);
% multi_comm_group_g2 = multilayer_community_detection_group(multi_comm_indivi_g2, ...
%     100, 'max', gamma);
% multi_comm_group_g3 = multilayer_community_detection_group(multi_comm_indivi_g3, ...
%     100, 'max', gamma);
% multi_comm_group_g4 = multilayer_community_detection_group(multi_comm_indivi_g4, ...
%     100, 'max', gamma);
% N_all_group_g1 = multi_comm_group_g1.multi_module_consensus;
% N_all_group_g2 = multi_comm_group_g2.multi_module_consensus;
% N_all_group_g3 = multi_comm_group_g3.multi_module_consensus;
% N_all_group_g4 = multi_comm_group_g4.multi_module_consensus;

%% 7e. Save Results
% output_dir = '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/modularity_var/';

cd '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/modularity_var/';

% filename = sprintf('S_indivi_g1_gamma.mat', gamma, log_omega); save(filename,'N_all_g1');
% % filename = sprintf('S_group_g1_gamma.mat', gamma, log_omega); save(filename,'N_all_group_g1');
% 
% filename = sprintf('S_indivi_g2_gamma.mat', gamma, log_omega); save(filename,'N_all_g2');
% % filename = sprintf('S_group_g2_gamma.mat', gamma, log_omega); save(filename,'N_all_group_g2');
% 
% filename = sprintf('S_indivi_g3_gamma.mat', gamma, log_omega); save(filename,'N_all_g3');
% % filename = sprintf('S_group_g3_gamma.mat', gamma, log_omega); save(filename,'N_all_group_g3');
% 
% filename = sprintf('S_indivi_g4_gamma.mat', gamma, omega); save(filename,'N_all_g4');
% % filename = sprintf('S_group_g4_gamma.mat', gamma, log_omega); save(filename,'N_all_group_g4');

%filename = sprintf('S_indivi_g_gamma.mat', gamma, omega); save(filename,'N_all_g');

% filename = sprintf('S_prefes.mat', gamma, log_omega); save(filename,'S_g1'); 
% filename = sprintf('S_prenfes.mat', gamma, log_omega); save(filename,'S_g2');
% filename = sprintf('S_postfes.mat', gamma, log_omega); save(filename,'S_g3');
% filename = sprintf('S_postnfes.mat', gamma, log_omega); save(filename,'S_g4');
% % filename = sprintf('S_g.mat', gamma, omega); save(filename,'S_g');

disp('Processing complete!');
