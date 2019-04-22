function groups = SpectralClustering(W, n, varargin) 
%SPECTRALCLUSTERING Spectral clustering.
%   GROUPS = SpectralClustering(W, N, varargin) clusters data into N groups
%   with affinity W.

% Input Arguments
% W                -- symmetric affinity matrix.
% n                -- number of groups.
% 'Start'          -- initial group for k-means.
%     'sample'(default):
%     the same as the k-means
% 'MaxIter'        -- maximum number of iterations for KMeans 
%     1000(default):
%     positive integer
% 'Replicates'     -- number of replications for KMeans
%     20(default):
%     positive integer
% 'Eig_Solver'     -- eig function of matlab
%     eig(default):
%     eigs

% Algorithm (Ncut)
% Ncut: min_{A_i} sum_i cut(A_i, bar(A)_i) / vol(A_i),
% reformulate as min_{A_i} \sum_i h_i' L h_i 
% s.t. h_ij = 1 / sqrt(Vol(A_j)), if v_i \in A_j
% relax as min_H tr(H' L H) s.t. H' D H = I
% solution is H = D^{-0.5} eig(L_sym).
% in this version, do another row normalization of H.

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

% check data
if ~issymmetric(W)
    error(['In ''' mfilename ''': affinity matrix is not symmetric'])
end
% define defaults options
% Set default 
vararg = {'Start', 'sample', ...
          'MaxIter', 1000, ...
          'Replicates', 20, ...
          'Eig_Solver', 'eig'};
% Overwrite by input
vararg = vararginParser(vararg, varargin);
% Generate variables
for pair = reshape(vararg, 2, []) % pair is {propName;propValue}
   eval([pair{1} '= pair{2};']);
end

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
% The computation is equivalent to:
% - compute the largest eigenvectors of D^{-1} W
% - normalize the rows of the resultant matrix
% - then apply kmeans to the rows.
if strcmpi(Eig_Solver, 'eig')
    [V, D] = eig( cnormalize(full(W), 1)' );
    [~, ord] = sort(real(diag(D)), 'descend');
    kerN = V(:, ord(1:n));
    clear V D;
elseif strcmpi(Eig_Solver, 'eigs')
    [kerN, ~] = eigs( cnormalize(W, 1)', n, 'LR' );
end
kerN = cnormalize_inplace(kerN')';
groups = kmeans(kerN, n, 'Start', Start, ...
                         'MaxIter', MaxIter, ...
                         'Replicates', Replicates, ...
                         'EmptyAction', 'singleton');
