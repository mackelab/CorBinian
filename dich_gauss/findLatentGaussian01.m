function [gamma,Lambda] = findLatentGaussian(mu,Sigma,acc)

% function [gamma,Lambda] = findLatentGaussian(mu,Sigma,acc)
% 	Compute mean gamma and covariance Lambda for the hidden Gaussian random vector Z generating the
% 	binary Bernoulli vector X with mean mu and covariances Sigma according to
%     X = 0 <==> Z < 0
%     X = 1 <==> Z > 0
%
% Usage: [gamma,Lambda] = findLatentGaussian([0.5 0.5],[0.25 0.1; 0.1 0.25])
%
% Code from the paper: 'Generating spike-trains with specified
% correlations', Macke et al., Neural Computation
%
%note: does NOT ensure that the matrix Lambda is positive definite, if oyu
%want to do this, need to use the higham-correction.


% parameter check
if any(mu < 0 | mu >= 1)
    error('mean has to be in (0,1)')
end

if nargin==2
    acc=1e-10;
end

% find mean and variances of the hidden variable via equ. 2.2
n = length(mu);
g = norminv(mu);
L = eye(n);

% find covariances by solving
% Sigma_ij - Psi(gamma_i,gamma_j,Lambda_ij) = 0 (equ. 2.2)
for i = 1:n
    for j = i+1:n
        cMin = -1;
        cMax = 1;
        
        % constant
        pn = prod(normcdf([g(i) g(j)]));
        
        %if covariance is 0, lets not do optimization and just write down
        %the answer:
        
        if abs(Sigma(i,j))<acc/2
            cMax=0;
            cMin=0;
        else %otherwise do optimization:
            
            % check whether DG distribution for covariance exists
            if (Sigma(i,j) - bivnor(-g(i),-g(j),-1) + pn) < -1e-3 || ...
                    (Sigma(i,j) - bivnor(-g(i),-g(j),1) + pn) > 1e-3
                error('A joint Bernoulli distribution with the given covariance matrix does not exist!')
                
            end
            
            % determine Lambda_ij iteratively by bisection (Psi is monotonic in rho)
            while cMax - cMin > acc
                cNew = (cMax + cMin) / 2;
                if Sigma(i,j) > bivnor(-g(i),-g(j),cNew) - pn
                    cMin = cNew;
                else
                    cMax = cNew;
                end
            end
        end
        L(i,j) = cMax;
        L(j,i) = cMax;
    end
end

gamma=g;
Lambda=L;

