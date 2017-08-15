function Efx = getEfxFromData(output, idxSubsamples, par, useSampleSizes )
% getEfxFromData
% Generate data means E[f(X)] over date X from (simulated) activity data 
% and indices specifying the subset of RGCs to consider. 

if nargin < 3
    useSampleSizes = 2:2:12; % n = 20, 40, 60, 80, 100, 120
end

% change format of activity data into something more intuitive
 n = size(output,3);         % get input spike raster format
 Nc = size(output,1);        % n = number of cells, Nc = trial length
 nTrials = size(output,2);   % nTrials = number of trials
 output = double(output);      % pyhton code produces data of type int64
 tmp = zeros(n, Nc, nTrials);
 for i = 1:nTrials
   tmp(:, :, i) = squeeze(output(:,i,:))';
 end
clear output;     % could pick a new name, but I got kind of to this one 
output = struct;  %
output.spikes = zeros(n, Nc*nTrials);
for i = 1:n
   output.spikes(i, :) = vec(squeeze(tmp(i,:,:))); % #cells-by-#data_points
end
clear Nc nTrials n tmp i

Efx = cell(length(idxSubsamples),1); % storage for data means
for i = useSampleSizes
  n = size(idxSubsamples{i},1); % current population size
  disp(['n = ', num2str(n)])
  Efx{i} = zeros(n*(n+3)/2+1, par.nDraws); % #features-by-#runs
  for j = 1:par.nDraws
    disp(['draw #', num2str(j), '/', num2str(par.nDraws)])
    tmp = zeros( n*(n+3)/2+1,1 ); % storage for running average
    for t = 1:floor(size(output.spikes,2)/100) % have to divide data into 
      tmp = tmp + full(mean( ...                % chunks to save memory
            setup_features_maxent(output.spikes(idxSubsamples{i}(:,j), ...
                                                (t-1)*100+(1:100))', ...
                                                 'ising_count_l_0'), 1)');
              
    end
    Efx{i}(:,j) = tmp/t; % we added t many means of data chunks of equal
                         % size, now re-normalize. Note that anything 
  end                    % between t*100 and size(output.spikes,2) got 
end                      % neglected - might be up to 99 date points!

end