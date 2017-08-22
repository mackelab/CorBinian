function [W, RGCcen] = genFilters(d,n,mode,pars)
% Input
%  - d: vector (2x1) or number (1x1) of stimuli dimensions. If vector, d(1)
%       and d(2) give extensions in y- and x-direction, respectively. If 
%       a number, the code assumes the stimuli to be rectangular with
%       dimensions [d,d] = size(stimulus);
%  - n: number of filters/RGCs to be generated
%  - mode: string defining the model for individual RFs and retinal tiling
%     - mode.RF:     supported values: 'center-surround DoG'
%     - mode.tiling: supported values: 'random', 'lattice',
%                                      'random-lattice'
%  - pars: additional, specific parameters for different models
%     - pars.thres:  filter values below threshold are set to zero
%     - pars.Sigma:  n-by-2 cell of covariance matrices for DoG filters, 
%                    or 1-by-2 cell of templates for all n DoG filters. 
%     - pars.hight:  n-by-2 matrix of respective hight of DoG components,
%                    or 1-by-2 matrix of templates for all n DoG filters.

if numel(d)==1 % can provide # of pixels both in X and Y or assume square
 if ~isinteger(sqrt(d))
   disp('Warning: Assuming square stimuli. Using rounded value');
 end
 d = [1;1]*ceil(sqrt(d));
end

if numel(pars.Sigma)==2 % if there is a single template for all DoG filters
  Sigma = cell(n,2);
  if iscell(pars.Sigma) % template is full covariance matrix
   for i = 1:n
     Sigma{n,1} = pars.Sigma{1};
     Sigma{n,2} = pars.Sigma{2};
   end
  else
   for i = 1:n          % template is just scale of unit matrix (iso-Gauss)
     Sigma{n,1} = diag(ones(2,1)*pars.Sigma{1});
     Sigma{n,2} = diag(ones(2,1)*pars.Sigma{2});
   end      
  end    
   pars.Sigma = Sigma; clear Sigma; 
end

if numel(pars.hight)==2 % if there is a single template for all DoG filters
  if size(pars.hight,1)==2
    pars.hight = pars.hight';
  end
  pars.hight = ones(n,1)*pars.hight; clear hight;
end

YX = [vec((ones(d(2),1) * (1:d(1)))'), vec((ones(d(1),1) * (1:d(2))))];

switch mode.tiling
    case 'random'
      RGCcen(1,:) = ((0.8*rand(1,n))+0.1)*d(1);
      RGCcen(2,:) = ((0.8*rand(1,n))+0.1)*d(2);
    case 'lattice'
      cy = round(linspace(1, d(1), ceil(sqrt(n))));
      cx = round(linspace(1, d(2), ceil(sqrt(n))));
      RGCcen(1,:) = vec(ones(length(cy),1) * cy);
      RGCcen(2,:) = vec((ones(length(cx),1) * cx)');
      RGCcen = RGCcen(:, 1:n);
end % end switch model.tiling

W = zeros(d(1)*d(2), n);
switch mode.RF 
    case 'center-surround DoG'
     for i = 1:n
       W(:,i) = pars.hight(i,1)*mvnpdf(YX,RGCcen(:,i)',pars.Sigma{i,1})...
              - pars.hight(i,2)*mvnpdf(YX,RGCcen(:,i)',pars.Sigma{i,2}); 
          
     end
     W = W';
     W(abs(W)<pars.thres) = 0;
end % end switch model.RF

end % end function


