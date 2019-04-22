function [X, Xnorm] = cnormalize_inplace(X, p)
%CNORMALIZE_INPLACE normalizes columns.
%   This is a inplace version of CNORMALIZE.

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

N = size( X, 2 );
if ~exist('p', 'var')
    p = 2;
end

if nargout > 1, Xnorm = zeros(1, N); end;

for iN = 1:N
    if p == Inf
        cnorm = max(abs(X(:, iN)), [], 1);
    else
        cnorm = sum(abs(X(:, iN)) .^p, 1) .^(1/p);
    end
    X(:, iN) = X(:, iN) / (cnorm + eps);
    
    if nargout > 1, Xnorm(iN) = cnorm; end;
end
