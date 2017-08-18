function [h,J,Z]=ReadResultDataFile(paramsresult,flipper)

fid = fopen(paramsresult, 'r');


% retrieve number of cells
try
ncells = fscanf(fid,'%d\n',1);
catch
    h=nan;
    J=nan;
    Z=nan;
    warning('results file not found, possibly the fitting did not complete successfully')
    return
end


% retrieve biases
h = fscanf(fid,'%g %g',[2 ncells]);



h = h(2,:)';

% retrieve couplings
Jdata = fscanf(fid,'%g %g %g',[3,inf])';
J = zeros(ncells,ncells);
J(sub2ind(size(J),Jdata(:,1)+1,Jdata(:,2)+1)) = Jdata(:,3);
J=2*J;
%J = J + J';

fclose(fid);

%
if nargin==2
    %do some flipping:
    A=diag(-2*flipper+1);
    J=A*J*A;
    %keyboard
    h=A*h-.5*(J+J')*flipper';
end

