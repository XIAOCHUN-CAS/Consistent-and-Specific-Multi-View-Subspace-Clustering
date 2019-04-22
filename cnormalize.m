function [Y, Xnorm] = cnormalize(X, p)
%CNORMALIZE normalizes columns.
% 
% [Y, XNORM] = cnormalize(X, P) normalizes columns of matrix X to unit
% ell_P norm, and returen the norm values to XNORM and data to Y.

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

if ~exist('p', 'var')
    p = 2;
end

if p == Inf
    Xnorm = max(abs(X), [], 1);
else
    Xnorm = sum(abs(X) .^p, 1) .^(1/p);
end
Y = bsxfun(@rdivide, X, Xnorm + eps);