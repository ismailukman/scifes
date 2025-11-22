function [B,twom] = multiord(A,gamma,omega)
%MULTIORD  returns multilayer Newman-Girvan modularity matrix for ordered layers, matrix version
% Works for directed or undirected networks
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
%   Input: A: Cell array of NxN adjacency matrices for each layer of an
%          ordered multilayer (directed or undirected) network
%          gamma: intralayer resolution parameter
%          omega: interlayer coupling strength
%
%   Output: B: [NxT]x[NxT] flattened modularity tensor for the
%           multilayer network with uniform ordinal coupling (T is
%           the number of layers of the network)
%           mm: normalisation constant
%
%   Example of usage: [B,mm]=multiord(A,gamma,omega);
%          [S,Q]= genlouvain(B); % see iterated_genlouvain.m and
%          postprocess_temporal_multilayer.m for how to improve output
%          multilayer partition
%          Q=Q/mm;
%          S=reshape(S,N,T);
%
%   [B,mm] = MULTIORD(A,GAMMA, OMEGA) with A a cell array of square
%   (symmetric or assymetric) matrices of equal size each representing a
%   directed or undirected network "layer" computes the Newman Girvan multilayer
%   modularity matrix using the quality function described in Mucha et al.
%   2010, with intralayer resolution parameter GAMMA, and with interlayer
%   coupling OMEGA connecting nearest-neighbor ordered layers.  The null
%   model used for the quality function is the Newman-Girvan null model
%   (see e.g. Bazzi et al. for other possible null models). Once the
%   mulilayer modularity matrix is computed, optimization can be performed
%   by the generalized Louvain code GENLOUVAIN or ITERATED_GENLOUVAIN. The
%   sparse output matrix B can be used with other heuristics, provided the
%   same mapping is used to go from the multilayer tensor to the multilayer
%   flattened matrix. That is, the node-layer tuple (i,s) is mapped to
%   i + (s-1)*N. [Note that we can define a mapping between a multilayer
%   partition S_m stored as an N by T matrix and the corresponding flattened
%   partition S stored as an NT by 1 vector. In particular S_m = reshape(S,N,T)
%   and S = S_m(:).]
%
%   See also
%       genlouvain heuristics:      GENLOUVAIN, ITERATED_GENLOUVAIN
%       multilayer wrappers:        MULTICAT, MULTICATF, MULTIORDF
%       other heuristics:           SPECTRAL23
%       Kernighan-Lin improvement:  KLNB
%
%   Notes:
%     The matrices in the cell array A are assumed to be square,
%     and of equal size.  These assumptions are not checked here.
%
%     This code assumes that the sparse quality/modularity matrix B will
%     fit in memory and proceeds to build that matrix.  For larger systems,
%     try MULTIORD_F for undirected layer networks and MULTIORDDIR_F
%     for directed layer networks.
%
%     This code serves as a template and can be modified for situations
%     with other wrinkles (e.g., different intralayer null models [see eg
%     Bazzi et al. 2016 for examples], different numbers of nodes from
%     layer-to-layer, or systems which are both multiplex and longitudinal).
%     That is, this code is only a starting point; it is by no means
%     exhaustive.
%
%     By using this code, the user implicitly acknowledges that the authors
%     accept no liability associated with that use.  (What are you doing
%     with it anyway that might cause there to be a potential liability?!?)
%
%   References:
%     Blondel, Vincent D., Jean-Loup Guillaume, Renaud Lambiotte, and
%     Etienne Lefebvre, "Fast unfolding of communities in large networks,"
%     Journal of Statistical Mechanics: Theory and Experiment, P10008
%     (2008).
%
%     Fortunato, Santo, "Community detection in graphs," Physics Reports
%     486, 75-174 (2010).
%
%     Good, Benjamin H., Yves-Alexandre de Montjoye, and Aaron Clauset,
%     "Performance of modularity maximization in practical contexts,"
%     Physical Review E 81, 046106 (2010).
%
%     Newman, Mark E. J. and Michelle Girvan. "Finding and Evaluating
%     Community Structure in Networks", Physical Review E 69, 026113 (2004).
%
%     Elizabeth A. Leicht and Mark E. J. Newman. "Community structure in
%     Directed Networks", Physical Review Letters 100, 118703 (2008).
%
%     Mucha, Peter J., Thomas Richardson, Kevin Macon, Mason A. Porter, and
%     Jukka-Pekka Onnela. "Community Structure in Time-Dependent,
%     Multiscale, and Multiplex Networks," Science 328, 876-878 (2010).
%
%     Bazzi, Marya, Mason A. Porter, Stacy Williams, Mark McDonald, Daniel
%     J. Fenn, and Sam D. Howison. "Community Detection in Temporal
%     Multilayer Networks, with an Application to Correlation Networks",
%     MMS: A SIAM Interdisciplinary Journal 14, 1-41 (2016).
%
%     Porter, M. A., J. P. Onnela, and P. J. Mucha, "Communities in
%     networks," Notices of the American Mathematical Society 56, 1082-1097
%     & 1164-1166 (2009).
%
%   Acknowledgments:
%     Thank you to Dani Bassett, Jesse Blocher, Bruce Rogers, and Simi Wang
%     for their collaborative help which led to significant cleaning up
%     of earlier versions of our multilayer community detection codes.

% Step 1: Set default parameters if not provided
if nargin<2
gamma=1; % Default resolution parameter = 1
end

if nargin<3
omega=1; % Default interlayer coupling = 1
end

% Step 2: Extract network dimensions
N=length(A{1});  % Number of nodes (brain regions) in each layer
T=length(A);     % Number of layers (time windows)


% Step 3: Handle scalar or vector gamma (allows different resolution per layer)
if length(gamma)==1
gamma=repmat(gamma,T,1);
end

% Step 4: Pre-allocate sparse modularity matrix B
% Size: (N×T) × (N×T) to hold all node-layer combinations
% Estimated non-zero elements: N²T (intralayer) + 2NT (interlayer)
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twom=0; % Initialize total edge weight (normalization constant)

% Step 5: Build intralayer modularity matrices for each layer
for s=1:T
    % Step 5a: Compute degree sequences for layer s
    kout=sum(A{s},1);  % Out-degree (row sum): [1 × N] vector
    kin=sum(A{s},2);   % In-degree (column sum): [N × 1] vector
    mm=sum(kout);      % Total edge weight in layer s
    twom=twom+mm;      % Accumulate total edges across all layers
    
    % Step 5b: Calculate node-layer indices for this layer
    % Maps (node i, layer s) to position i + (s-1)*N in flattened matrix
    indx=[1:N]+(s-1)*N;
    
    % Step 5c: Build Newman-Girvan modularity matrix for layer s
    % Formula: B_ij = (A_ij + A_ji)/2 - gamma * (k_i*k_j + k_j*k_i)/(2m)
    % This penalizes edges expected under configuration model
    B(indx,indx)=(A{s}+A{s}')/2-gamma(s)/2.*((kin*kout+kout'*kin')/mm);
end

% Step 6: Add interlayer coupling (NEAREST-NEIGHBOR ONLY)
% Creates diagonal offsets connecting ONLY ADJACENT layers
% For T=3: connects layer 1-2, 2-3 (but NOT 1-3)
% Offset -N: connects layer s to layer s-1
% Offset +N: connects layer s to layer s+1
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);

% Step 7: Update normalization constant with interlayer coupling
% Only adjacent layers connected: (T-1) pairs of connections
% Each pair has 2N connections (N nodes × 2 directions)
% Total interlayer connections: 2*N*(T-1)
twom=twom+2*N*(T-1)*omega;
