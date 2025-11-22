function [B,twom] = modularity(A,gamma)
%MODULARITY returns monolayer Newman-Girvan modularity matrix for network given by adjacency matrix A, matrix version
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
%Works for directed and undirected networks
%
%   Input: A:  NxN adjacency matrices of a directed or undirected network
%          gamma: resolution parameter
%
%   Output: B: function handle where B(i) returns the ith column of
%          [N]x[N] modularity matrix of the monolayer network
%           with adjacency matrix A
%           twom: normalisation constant
%
%   Example of usage: [B,twom]=modularity(A,gamma);
%          [S,Q]= genlouvain(B);
%          Q=Q/twom;
%   Notes:
%     The matrix A is assumed to be square. This assumption is not checked
%     here.
%
%     This code assumes that the sparse quality/modularity matrix B will
%     fit in memory and proceeds to build that matrix.  For larger systems,
%     try MODULARITY_F for undirected networks and MODULARITYDIR_F for directed
%     networks.
%
%     This code serves as a template and can be modified for situations
%     with other wrinkles (e.g., different null models).
%
%     By using this code, the user implicitly acknowledges that the authors
%     accept no liability associated with that use.  (What are you doing
%     with it anyway that might cause there to be a potential liability?!?)
%
%   References:
%     Newman, Mark E. J. and Michelle Girvan. "Finding and Evaluating
%     Community Structure in Networks", Physical Review E 69, 026113 (2004).

fprintf('modularity from genlouvain... \n');

% Step 1: Set default resolution parameter if not provided
if nargin<2||isempty(gamma)
	gamma=1; % Default resolution parameter = 1
end

% Step 2: Compute degree sequences from adjacency matrix 
k=sum(A,2);  % In-degree: sum across columns, gives [N × 1] vector
             % For undirected networks, this equals node degree

d=sum(A,1);  % Out-degree: sum across rows, gives [1 × N] vector
             % For undirected networks, this also equals node degree

% Step 3: Calculate total edge weight (normalization constant)
twom=sum(k);  % Sum of all degrees = 2m (total edges × 2 for undirected)
              % For directed networks, this is total edge weight

% Step 4: Build Newman-Girvan modularity matrix
% Formula: B_ij = A_ij - gamma * (k_i * d_j) / (2m)
% Components:
%   - (A + A')/2: Symmetrize adjacency matrix for undirected networks
%   - k*d: Expected edges from configuration model (outer product [N×1] × [1×N] = [N×N])
%   - d'*k': Transpose to ensure symmetry
%   - gamma: Resolution parameter (higher = more/smaller communities)
%   - twom: Normalization by total edge weight
B=full((A+A')/2-gamma/2*(k*d+d'*k')/twom);

end
