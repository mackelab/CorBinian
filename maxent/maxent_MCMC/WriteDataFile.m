function A=WriteDataFile(datafile,mu,cov)
% writes mu and cov data into datafile
A=[];
mu=mu(:)';

fout = fopen(datafile, 'wt');
ncells = length(mu);
[r c] = size(cov);
if (r ~= c || ncells ~= r)
    disp('Error: mu and cov sizes do not match!');
end

fprintf(fout, '%d\n',ncells);
mu = [0:(ncells-1) ; mu]';
fprintf(fout, '%d %15.14f\n', mu');

subs = nchoosek(1:ncells,2);
inds = sub2ind(size(cov), subs(:,1), subs(:,2));
fprintf(fout, '%d %d %15.14f\n', [subs-1, cov(inds)]');
fclose(fout);
