function extractKeyStats(data, Nc, figNumber)

if nargin < 3 || isempty(figNumber)
    figNumber = 42^2;
end
if size(data.spikes,2) < size(data.spikes,1)
  data.spikes = data.spikes';
end

if mod(size(data.spikes,2),Nc)~=0
    error('data not consistent with N repetitions of a stimulus of length Nc')
end

n = size(data.spikes,1);
N = size(data.spikes,2);
data.spikes = reshape(data.spikes, [n, Nc, N/Nc]);

cRange = -0.1:0.01:0.6;
hFull = histc(data.spkCorrs(logical(triu(ones(n),1))), cRange);
hFull = hFull/sum(hFull);
disp('Extracting signal correlations (from averages)')
data.sigCorrs = corr(squeeze(mean(data.spikes,3))');
hSig = histc(data.sigCorrs(logical(triu(ones(n),1))), cRange);
hSig = hSig/sum(hSig);
disp('Extracting signal correlations (from shuffling)')
tmp = zeros(size(data.spikes));
for i = 1:n
    idx = randperm(N/Nc);
    tmp(i,:,:) = data.spikes(i, :, idx);
end
tmp = squeeze(reshape(tmp, [n, N, 1]));
data.sigShuffCorrs = corr(tmp');
hSigShuff = histc(data.sigShuffCorrs(logical(triu(ones(n),1))), cRange);
hSigShuff = hSigShuff/sum(hSigShuff);
hNoi = histc(data.spkCorrs(logical(triu(ones(n),1)))-data.sigShuffCorrs(logical(triu(ones(n),1))), cRange);
hNoi = hNoi/sum(hNoi);
M = max([hFull(:);hSig(:);hSigShuff(:);hNoi(:)]);

figure(figNumber);
disp('Extracting full correlations')
subplot(1,3,1)
bar(cRange, hFull);
title(['full correlations corr\_full, ', num2str(n), ...
       ' neurons, ', num2str(N/Nc), ' trials'])
xlabel('correlation coeff'); ylabel('emp. distribution of corr. coeff.')
box off; set(gca, 'TickDir', 'out')
axis([min(cRange)-(cRange(2)-cRange(1))/2, max(cRange)+(cRange(2)-cRange(1))/2, -0.05*M, 1.05*M]);
text(0.5*max(cRange), 0.3*M, ['std[x]: ',num2str(round(sqrt(var(data.spkCorrs(logical(triu(ones(n),1)))))*1000)/1000)]);
text(0.5*max(cRange), 0.6*M, ['E[|x|]: ',num2str(round(mean(abs(data.spkCorrs(logical(triu(ones(n),1)))))*1000)/1000)]);

% subplot(2,2,2)
% bar(cRange, hSig);
% title(['correlations in PSTH, corr_psth, ', num2str(n), ...
%        ' neurons, ', num2str(N/Nc), ' trials'])
% xlabel('correlation coeff'); ylabel('emp. distribution of corr. coeff.')
% box off; set(gca, 'TickDir', 'out')
% axis([min(cRange)-(cRange(2)-cRange(1))/2, max(cRange)+(cRange(2)-cRange(1))/2, -0.05*M, 1.05*M]);
% text(0.5*max(cRange), 0.3*M, ['std[x]: ',num2str(round(sqrt(var(data.sigCorrs(logical(triu(ones(n),1)))))*1000)/1000)]);
% text(0.5*max(cRange), 0.6*M, ['E[|x|]: ',num2str(round(mean(abs(data.sigCorrs(logical(triu(ones(n),1)))))*1000)/1000)]);

subplot(1,3,2)
bar(cRange, hSigShuff);
title(['signal correlations corr\_sig, ', num2str(n), ...
       ' neurons, ', num2str(N/Nc), ' shuffled trials'])
xlabel('correlation coeff'); ylabel('emp. distribution of corr. coeff.')
box off; set(gca, 'TickDir', 'out')
axis([min(cRange)-(cRange(2)-cRange(1))/2, max(cRange)+(cRange(2)-cRange(1))/2, -0.05*M, 1.05*M]);
text(0.5*max(cRange), 0.3*M, ['std[x]: ',num2str(round(sqrt(var(data.sigShuffCorrs(logical(triu(ones(n),1)))))*1000)/1000)]);
text(0.5*max(cRange), 0.6*M, ['E[|x|]: ',num2str(round(mean(abs(data.sigShuffCorrs(logical(triu(ones(n),1)))))*1000)/1000)]);

subplot(1,3,3)
bar(cRange, hNoi);
title(['noise correlations corr\_full - corr\_sig, ', num2str(n), ...
       ' neurons'])
xlabel('correlation coeff'); ylabel('emp. distribution of corr. coeff.')
box off; set(gca, 'TickDir', 'out')
axis([min(cRange)-(cRange(2)-cRange(1))/2, max(cRange)+(cRange(2)-cRange(1))/2, -0.05*M, 1.05*M]);
text(0.5*max(cRange), 0.3*M, ['std[x]: ',num2str(round(sqrt(var(data.spkCorrs(logical(triu(ones(n),1)))-data.sigShuffCorrs(logical(triu(ones(n),1)))))*1000)/1000)]);
text(0.5*max(cRange), 0.6*M, ['E[|x|]: ',num2str(round(mean(abs(data.spkCorrs(logical(triu(ones(n),1)))-data.sigShuffCorrs(logical(triu(ones(n),1)))))*1000)/1000)]);


figure(figNumber+1)
subplot(1,3,1)
imagesc(data.spkCorrs); colorbar
title('Full correlations')
% subplot(2,2,2)
% imagesc(data.sigCorrs); colorbar
% title('(averaged) signal correlations')
subplot(1,3,2)
imagesc(data.sigShuffCorrs); colorbar
title('(shuffled) signal correlations')
subplot(1,3,3)
imagesc(data.spkCorrs-data.sigShuffCorrs); colorbar
title('full - (shuffled) signal correlations')

end
