clear
clc

gamma = 1.21; %1.2 gamma = 0.5:0.1:1.5
log_omega = -1; %-1 log_omega = 0:-1:-4
omega = power(10,log_omega);
% omega = 0.1;

corr_g1 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_fes.mat'))), [2 3 1]);
corr_g2 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_pre_nfes.mat'))), [2 3 1]);
corr_g3 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/fc/corr_post_fes.mat'))), [2 3 1]);

% corr_g1 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_GraphTheory/sci_data/SCI/fc/corr_sci_c.mat'))), [2 3 1]);
% corr_g2 = permute(cell2mat(struct2cell(load('/Users/ismaila/Documents/C-Codes/SCI_GraphTheory/sci_data/SCI/fc/corr_sci_t.mat'))), [2 3 1]);


% model settings 
A_g1 = squeeze(num2cell(corr_g1,[1 2])); % adjacency/connectivity matrix
A_g2 = squeeze(num2cell(corr_g2,[1 2]));    
A_g3 = squeeze(num2cell(corr_g3,[1 2]));
N = length(A_g1{1}); % number of nodes (same for both groups)
T_g1 = size(corr_g1,3); % number of layers (subjects)
T_g2 = size(corr_g2,3);
T_g3 = size(corr_g3,3);

% multi-layer modularity (Q) and partitions
% Group 1
[B_g1, twom_g1] = multicat(A_g1,gamma,omega);
% [B_g1, twom_g1] = multicatdir_f(A_g1,gamma,omega);
PP_g1 = @(S_g1)postprocess_categorical_multilayer(S_g1,T_g1);
% [S_g1,Q_g1] = iterated_genlouvain(B_g1,10000,0,1,'moverandw',[], PP_g1); % 4th entry(randord): 0[move] or 1[moverand] | 5th: move, moverand, or moverandw
[S_g1, Q_g1] = iterated_genlouvain(B_g1, 10000, 0, 1, 'moverandw', [], PP_g1);
Q_g1 = Q_g1/twom_g1; % Community Modularity
% fprintf('Modularity Group 1: %d\n', Q_g1);
S_g1 = reshape(S_g1,N,T_g1);
C_g1 = mode(S_g1,2); % consensus
K_g1 = max(S_g1,[],'all'); % number of communities

% Group 2
[B_g2,twom_g2] = multicat(A_g2,gamma,omega);
PP_g2 = @(S_g2)postprocess_categorical_multilayer(S_g2,T_g2);
% [S_g2,Q_g2] = iterated_genlouvain(B_g2,10000,0,1,'moverandw',[], PP_g2); % 4th entry(randord): 0[move] or 1[moverand] | 5th: move, moverand, or moverandw
[S_g2, Q_g2] = iterated_genlouvain(B_g2, 10000, 0, 1, 'moverandw', [], PP_g2);
Q_g2 = Q_g2/twom_g2; % Community Modularity
% fprintf('Modularity Group 2: %d', Q_g2);
S_g2 = reshape(S_g2, N, T_g2);
C_g2 = mode(S_g2, 2); % consensus
K_g2 = max(S_g2,[],'all'); % number of communities

% Group 3
[B_g3,twom_g3] = multicat(A_g3,gamma,omega);
PP_g3 = @(S_g3)postprocess_categorical_multilayer(S_g3,T_g3);
% [S_g2,Q_g2] = iterated_genlouvain(B_g2,10000,0,1,'moverandw',[], PP_g2); 
[S_g3, Q_g3] = iterated_genlouvain(B_g3, 10000, 0, 1, 'moverandw', [], PP_g3);
Q_g3 = Q_g3/twom_g3; % Community Modularity
% fprintf('Modularity Group 2: %d', Q_g2);
S_g3 = reshape(S_g3, N, T_g3);
C_g3 = mode(S_g3, 2); % consensus
K_g3 = max(S_g3,[],'all'); % number of communities

%       'moverandw': move the node under consideration to a community chosen
%           at random from all moves that increase the qualilty where the
%           probability of choosing a particular move is proportional to
%           its increase in the quality function
% https://github.com/GenLouvain/GenLouvain/blob/master/iterated_genlouvain.m
% https://brainconn.readthedocs.io/en/latest/generated/brainconn.modularity.modularity_louvain_und.html

cd '/Users/ismaila/Documents/C-Codes/SCI_FES_GraphAnalysis/sci_data/SCI/modularity_var/';

filename = sprintf('S_prefes_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g1');
filename = sprintf('S_prenfes_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g2');
filename = sprintf('S_postfes_%.1f,%.1f.mat', gamma, log_omega); save(filename,'S_g3');


