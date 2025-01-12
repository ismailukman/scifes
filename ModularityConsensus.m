clear
clc

% corr_g1 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_GraphTheory/sci_data/SCI/fc/corr_hc.mat'))), [2 3 1]);
% corr_g2 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_GraphTheory/sci_data/SCI/fc/corr_sci.mat'))), [2 3 1]);

corr_g1 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_fes.mat'))), [2 3 1]);
corr_g2 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_nfes.mat'))), [2 3 1]);
corr_g3 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_fes.mat'))), [2 3 1]);

tic
n = 50;

% m_g1 = zeros(1, n);
% m_g2 = zeros(1, n);
m_g1 = cell(1, n);
m_g2 = cell(1, n);
m_g3 = cell(1, n);

k_g1 = cell(1, n);
k_g2 = cell(1, n);
k_g3 = cell(1, n);

gamma = 1.21; %1.2 gamma = 0.5:0.1:1.5
log_omega = -1; %-1 log_omega = 0:-1:-4
omega = power(10,log_omega);  

parfor i = 1:n
    % Use temporary variables for x and y inside parfor
    
    A_g1 = squeeze(num2cell(corr_g1, [1 2]));
    A_g2 = squeeze(num2cell(corr_g2, [1 2]));
    A_g3 = squeeze(num2cell(corr_g3, [1 2]));
    N = length(A_g1{1});
    T_g1 = size(corr_g1, 3);
    T_g2 = size(corr_g2, 3);
    T_g3 = size(corr_g3, 3);

    [B_g1, twom_g1] = multicat(A_g1,gamma,omega);
    PP_g1 = @(S_g1)postprocess_categorical_multilayer(S_g1, T_g1);
    [S_g1, Q_g1] = iterated_genlouvain(B_g1, 10000, 0, 1, 'moverandw', [], PP_g1);
    Q_g1 = Q_g1 / twom_g1;
    S_g1 = reshape(S_g1, N, T_g1);
    C_g1 = mode(S_g1, 2);
    K_g1 = max(S_g1, [], 'all');

    [B_g2, twom_g2] = multicat(A_g2, gamma, omega);
    PP_g2 = @(S_g2)postprocess_categorical_multilayer(S_g2, T_g2);
    [S_g2, Q_g2] = iterated_genlouvain(B_g2, 10000, 0, 1, 'moverandw', [], PP_g2);
    Q_g2 = Q_g2 / twom_g2;
    S_g2 = reshape(S_g2, N, T_g2);
    C_g2 = mode(S_g2, 2);
    K_g2 = max(S_g2, [], 'all');

    [B_g3, twom_g3] = multicat(A_g3, gamma, omega);
    PP_g3 = @(S_g3)postprocess_categorical_multilayer(S_g3, T_g3);
    [S_g3, Q_g3] = iterated_genlouvain(B_g3, 10000, 0, 1, 'moverandw', [], PP_g3);
    Q_g3 = Q_g3 / twom_g3;
    S_g3 = reshape(S_g3, N, T_g3);
    C_g3 = mode(S_g3, 2);
    K_g3 = max(S_g3, [], 'all');
    
    m_g1{i} = S_g1;
    m_g2{i} = S_g2;
    m_g3{i} = S_g3;

    k_g1{i} = K_g1;
    k_g2{i} = K_g2;
    k_g3{i} = K_g3;

    fprintf('%d. ',i);        
end

S_g1_temp = cat(3, m_g1{:});
S_g2_temp = cat(3, m_g2{:});
S_g3_temp = cat(3, m_g3{:});

k_g1_mode = mode(cell2mat(k_g1));
k_g2_mode = mode(cell2mat(k_g2));
k_g3_mode = mode(cell2mat(k_g2));
% 
% k_g1_mode = mode(k_g1);
% k_g2_mode = mode(k_g2);

S_g1 = S_g1_temp(:,:,5);
S_g2 = S_g2_temp(:,:,5);
S_g3 = S_g3_temp(:,:,5);
% find the mode across 3rd axis
S_g1_mode = mode(S_g1_temp, 3);
S_g2_mode = mode(S_g2_temp, 3);
S_g3_mode = mode(S_g3_temp, 3);

% find the average across 3rd axis
% S_g1_avg = mean(S_g1_temp, 3);
% S_g2_avg = mean(S_g2_temp, 3);

% S_g1_avg = round(mean(S_g1_temp, 3));
% S_g2_avg = round(mean(S_g2_temp, 3));

toc

cd '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/modularity_var/';

% filename = sprintf('S_hc_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g1');
% filename = sprintf('S_sci_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g2');
% filename = sprintf('Sa_hc_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g1_avg');
% filename = sprintf('Sa_sci_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g2_avg');

% filename = sprintf('Sm_hc_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g1_mode');
% filename = sprintf('Sm_sci_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g2_mode');
% filename = sprintf('Sm_scic_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g1_mode');
% filename = sprintf('Sm_scit_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g2_mode');

% filename = sprintf('Sm_hc_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g1_temp');
% filename = sprintf('Sm_sci_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g1_temp');

% filename = sprintf('var_s1.mat'); save(filename,'S_g1_temp');
% filename = sprintf('var_s2.mat'); save(filename,'S_g2_temp');
% filename = sprintf('var_s1_mode.mat'); save(filename,'S_g1_mode');
% filename = sprintf('var_s2_mode.mat'); save(filename,'S_g2_mode');

filename = sprintf('var_s1_g1.mat'); save(filename,'S_g1_temp');
filename = sprintf('var_s2_g1.mat'); save(filename,'S_g2_temp');
filename = sprintf('var_s3_g1.mat'); save(filename,'S_g3_temp');

filename = sprintf('var_s1_mode.mat'); save(filename,'S_g1_mode');
filename = sprintf('var_s2_mode.mat'); save(filename,'S_g2_mode');
filename = sprintf('var_s3_mode.mat'); save(filename,'S_g3_mode');

fprintf('Done!\n');


% % Finding the probability of same community assigning of node i, j
% % Initialize variables
% totalIterations = 100;
% nodeCount = 200;
% sameCommunityCount = 0;
% 
% % Loop over iterations
% for iter = 1:totalIterations
%     % Generate random community assignments for each node
%     communityAssignments = detectCommunities(nodeCount);
% 
%     % Check if nodes i and j are assigned to the same community
%     % (Assuming i and j are specific nodes of interest)
%     iCommunity = communityAssignments(i);
%     jCommunity = communityAssignments(j);
% 
%     if iCommunity == jCommunity
%         sameCommunityCount = sameCommunityCount + 1;
%     end
% end
% 
% % Calculate probability
% probability = sameCommunityCount / totalIterations;
% disp(['Probability of nodes i and j being in the same community: ' num2str(probability)]);
