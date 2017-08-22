data
clear all
close all
clc

% Settings concerning the overall retina simulation
%--------------------------------------------------------------------------
d = [50, 50]; % dimension of visual field, in pixels
n = 40;       % number of RGCs to be simulated
nComp = 40;   % number of RGCs to be included in a computational run (sub-sample)
idxC = randsample(n,nComp); % index of nComp cells chosen for a run
pRet.Ce=0.25*eye(d(1)*d(2)); % covariance matrix for Gaussian-induced noise correlations 
                             % currently white. Should reflect spatial
                             % correlations eventually. 
pRet.magnitude =  1;      % parameters governing the nonlinearity mapping
pRet.gain      =  1;      % from linear filter responses to RGC output
pRet.offset    = -20;     % spiking probability (sigmoidal)
mode.tiling = 'random';           % arrangement of RGC centres 
mode.RF = 'center-surround DoG';  % layout of RGC receptive fields

% Settings concerning input images to be fed into the model retina
%--------------------------------------------------------------------------
N = 50000; % number of image frames to be generated
Nc = 1000;  % chunks of images to be generated at a time (think of video)
alpha =  2; % 0 'should' be white, 1 pink and 2 brown noise 
xMag  =  1000; % magnitude of Gaussian noise to generate images from

% Parameters concerning the specific layout of the RGC receptive fields
%--------------------------------------------------------------------------
pars.thres = 0.00001;% threshold below which the RFs will be truncated to 0
Sigma = { [4,   0 ;      % template for central part of DoG filters
            0,  4 ], ... % (peak of Mexican heat)
          [9,   0 ;      % template for outer part of DoG filters
            0,  9 ]};    % (surround of Mexican hat)
SigmaDFs = 100000*[1,1];  % degrees of freeom for DoG component covariances
ON_OFF = 2*randi(2, [n,1]) - 3; % per cell whether it's an ON or an OFF RGC
hight = [1, 1]; % template for hights of central and outer DoG components
hightSTDs = [0.01, 0.01]; % standard deviates for hights of DoG components 
idxON  = find(ON_OFF > 0); % quick lists for finding
idxOFF = find(ON_OFF < 0); % ON and OFF cells
pars.hight = zeros(n,2); % hights of DoG component Gaussian bumps
pars.hight(:,1) = abs(normrnd(hight(1),hightSTDs(1)^2,[n,1])) .*ON_OFF;
pars.hight(:,2) = abs(normrnd(hight(2),hightSTDs(2)^2,[n,1])) .*ON_OFF;
for i = 1:n 
  pars.Sigma{i,1} = wishrnd(Sigma{1}, SigmaDFs(1))/SigmaDFs(1);
  pars.Sigma{i,2} = wishrnd(Sigma{2}, SigmaDFs(2))/SigmaDFs(2);
end

%%------------------------------------------------------------------------
% (2) Generate input images and do the retina simulation
%%-------------------------------------------------------------------------
disp('Generating RGC filters')
[W, RGCcen] = genFilters(d,n,mode,pars); 
W=sparse(W);
%--------------------------------------------------------------------------
disp('Visualizing RGC filters')
figure(1); % 3D images showing individual filters with different colors
subplot(1,2,1)
for i = 1:10; % ON cells 
 xy = reshape(W(idxON(i),:), d(1), d(2)); 
 xy(xy==0) = NaN;                  % don't plot these pixels, just clutters
 surf(xy, ones(d)*i/length(idxON)), hold on, 
end;
hold off, title('ON RGCs')
subplot(1,2,2)
for i = 1:length(idxOFF); % OFF cells
 xy = reshape(W(idxOFF(i),:), d(1), d(2)); 
 xy(xy==0) = NaN;                  % don't plot these pixels, just clutters
 surf(xy, ones(d)*i/length(idxOFF)), hold on, 
end;
hold off, title('OFF RGCs')
figure(2); % 2D image showing sum of filters, i.e. overall retina imprint
subplot(1,2,1), xy = reshape(sum(W(idxON,:),1), d(1), d(2));
imagesc(xy), title('ON RGCs (sum of all filters, check for full coverage mostly')
subplot(1,2,2), xy = reshape(sum(W(idxOFF,:),1), d(1), d(2));
imagesc(xy), title('OFF RGCs')
clear xy
% subplot(1,2,2); alternatively show centers of RGC mosaic
% plot(RGCcen(1,:), RGCcen(2,:), 'r.'); hold on;
% plot(RGCcen(1,:), RGCcen(2,:), 'ro', 'markerSize', 50);

disp('Generating input images and simulating RGC output')
out.spikes = zeros(n,N);
x = xMag * spatialPattern([d(1),d(2),Nc], -alpha);
for i = 1:floor(N/Nc)
  disp([' - chunk ', num2str(i), ' out of ', num2str(floor(N/Nc))])
  tmp = retSim(x,W,pRet);
  out.spikes(:,(i-1)*Nc+1:i*Nc) = tmp.spikes;
end
out.spkCorrs = full(corr(out.spikes'));
out.spkCov   = full(cov(out.spikes'));
out.RFoverlap = full(W*W');
%--------------------------------------------------------------------------
disp('Visualizing neural activity')
figure(3); 
subplot(2,3,1:3), imagesc(out.spikes), title('spike trains')
xlabel('time t'), ylabel('neuron n')
set(gca, 'tickdir', 'out'), box off;
subplot(234), imagesc(out.spkCorrs-diag(diag(out.spkCorrs)))
title('output spike correlations')
xlabel('neuron n'), ylabel('neuron n')
set(gca, 'tickdir', 'out'), box off;
subplot(235), imagesc(out.RFoverlap-diag(diag(out.RFoverlap)))
title('receptive field "overlaps" w(:,i)^T * w(:,j)')
xlabel('neuron n'), ylabel('neuron n')
set(gca, 'tickdir', 'out'), box off;
idxM = (logical(triu(ones(n,n),1))); 
OnOff = logical((ON_OFF*ON_OFF'+1)/2);
subplot(236), 
plot(out.spkCorrs(idxM&OnOff), out.RFoverlap(idxM&OnOff),'r.');
hold on
plot(out.spkCorrs(idxM&~OnOff), out.RFoverlap(idxM&~OnOff),'g.');
title('spike correlations = f(RF overlap)?')
xlabel('receptive field overlap')
ylabel('output spike correlations')
legend('ON-ON or OFF-OFF pairs', 'ON-OFF or OFF-ON pairs', 'location', 'southeast')
set(gca, 'tickdir', 'out'), box off;


figure(4); 
subplot(1,2,1);
h = histc(out.spkCorrs(logical(triu(ones(nComp),1))), -0.2:0.01:0.6); 
xlabel('corr. coeff.')
ylabel('relative frequency')
bar(-0.2:0.01:0.6, h/sum(h))
title('Distribution of correlation coefficients')
subplot(1,2,2);
plot(mean(out.spikes,2));
title('Distribution of firing rates (per bin, i.e. think of *20 for 50ms bins')

% slightly more elaborate than figure 4....
extractKeyStats(out, Nc, 5);