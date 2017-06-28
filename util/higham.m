function [A val it] = higham(B,acc,maxit)
% function A = higham(B)
%   Converts the matrix B which is assumed to be symetric to a valid
%   correlation matrix A, i.e. A_ii = 1, A positive semidefinite
%   using NJ Higham's iterative projection algorithm. The result is
%   the unique matrix X minimizing ||A-X||, i.e. the closest matrix in the
%   Frobenius norm sense.
%
%   Reference: NJ Higham, Computing the nearest correlation matrix - a
%   problem from finance, IMA Journal of Numerical Analysis, 2002
%
% Code from the paper: 'Generating spike-trains with specified
% correlations', Macke et al., Neural Computation
%
% www.kyb.mpg.de/bethgegroup/code/efficientsampling

% set defaults
if nargin < 3
    maxit = 1000;
end
if nargin < 2
    acc = 1e-6;
end

% initialization
it = 0; DS = 0;
Yold = B; Xold = B;
val = inf;

% iterate while changing or until iteration limit is reached
while val > acc && it < maxit
    % apply dykstra correction
    R = Yold - DS;
    
    % project onto S (positive definite matrices)
    Xnew = P_S(R);
    DS = Xnew - R;
    
    % project onto U (matrices with unit diagonal)
    Ynew = P_U(Xnew);
    
    % compute change
    nx = norm(Xnew-Xold,'inf')/norm(Xnew,'inf');
    ny = norm(Ynew-Yold,'inf')/norm(Ynew,'inf');
    nxy = norm(Ynew-Xnew,'inf')/norm(Ynew,'inf');
    val = max([nx ny nxy]);
    
    % update variables
    Xold = Xnew;
    Yold = Ynew;
    
    it = it + 1;
end
A = Ynew;

% throw warning if not converged
if it == maxit
    warning(myErrorMsg('Iteration limit reached without convergence.', ...
        [],'higham','higham.m',44))
    A = higham(A,1e-10,1e5);
end

%ensure that matrix is positive definite in the end:

mineig=min(eig(A));
if mineig<0
    warning('Somehow, correlation matrix was not pd, applying correction')
    [a,b]=eig(A);
    b=diag(max(diag(b),1e-8));
    A=a*b*a';
    A=cov_2_corr(A);
    A=A/2+A/2';
end


end


function X = P_U(Y)
t = diag(Y-eye(size(Y)));
X = Y - diag(t);
end

function X = P_S(Y)
[V,D] = eig(Y);
D(D<0) = 0;
X = V*D*V';
end
