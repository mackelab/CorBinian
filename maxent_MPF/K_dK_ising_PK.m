function [K, dK] = K_dK_ising_PK( lambda, data )
% Objective function for minimum probability flow model fitting for Ising
% models with additional terms ('K-pairwise maximum entropy models')

% Original author: Jascha Sohl-Dickstein (2012)
% Web: http://redwood.berkeley.edu/wiki/Jascha_Sohl-Dickstein

% Inputs: 
% lambda = [ J(:),L(:) ]: Concatentation of extended Ising+V(K) parameters
% data.x: d-by-n matrix, columns are binary data vectors
% data.counts = sum(data.x,1);
% data.mask: d-by-n boolean matrix, gives which of the bit-flipped
%            data vectors are already present in the raw data.

% Outputs:
% K:  K(lambda), the evaluated objective function of MPF
% dK: gradient, first n^2 elements correspond to J, rest to L

    % Precomputations
    %------------------------------------------------------------
    [d, n] = size( data.x );

    J = reshape(lambda(1:d^2), d, d);
    L = lambda(end-d:end);                
      %  assumes size(L) = [d+1,1], i.e. also a feature for K = 0!!!

    
    % Compute objective function
    %------------------------------------------------------------
    diagJ = diag(J); 
    dxn   = 2*data.x-1; % in n-th row flips n-th entries of each x 
    % Kfull is a [d, n] matrix containing the contribution to the
    % objective function from flipping each bit in the rows
    k   = data.counts(ones(d,1),:) + 1; % activity counts (+1 for indexing)
    Kfull = exp(dxn .* (J*data.x) +  ...
                (1/2) * ( -diagJ(:,ones(1,n)) + L(k) - L(k-dxn) ) ); 
    
    if isfield(data, 'doubleMask') % correction if samples differ by one bit
      Kfull = Kfull .* data.doubleMask;
    end
    
    K = sum(Kfull(:)); 
    
    
    % Compute derivatives dK/dJ of the standard Ising parameters 
    %------------------------------------------------------------
    dJ = (Kfull.*dxn) * data.x' - (1/2)*diag( sum(Kfull, 2) );
    dJ = (dJ + dJ')/2;
 
    
    % Compute derivatives dJ/dL of the activity count extension
    %------------------------------------------------------------
    indCount = logical(full(sparse(1:n,data.counts+1,1,n,d+1)));
    indCountMinu1 = zeros(n, d+1);                % Using linear algebra
    indCountMinu1(:,2:end) = indCount(:,1:end-1); % notation for simplicity
    indCountPlus1 = zeros(n, d+1);                % (and Matlab's somewhat
    indCountPlus1(:,1:end-1) = indCount(:,2:end); % fast). Could improve...
    
    dL =  (1/2) * ( sum(         Kfull,         1) * indCount ...
                  - sum(   data.x   .*   Kfull, 1) * indCountPlus1 ...
                  - sum( (1-data.x) .*   Kfull, 1) * indCountMinu1 );
    
    % Assemble output
    %------------------------------------------------------------
    K  = K  / n;
    dK = [dJ(:); dL(:)];
    dK = dK / n;
    
