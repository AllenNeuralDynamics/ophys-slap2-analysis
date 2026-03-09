function [skIm, P] = localizeSources_vIM(IM, vIM, params, doPlot)
%inputs:
%IM:        3D recording, X x Y x Time
%aData:     alignment metadata
nTimePoints = size(IM,3);
tau = params.tau_s.*params.alignHz; %time constant in frames
params.tau_frames = tau;
sigma = params.sigma_px; %space constant in pixels
baselineWindow = ceil(params.baselineWindow_Glu_s.*params.alignHz);
denoiseWindow = ceil(params.denoiseWindow_s.*params.alignHz);
nans = isnan(IM);

sz = size(IM);
valid = mean(nans,3)<params.nanThresh; %a pixel must be imaged at least (1-nanThresh) of the time to be included
if ~any(valid)
    warning('Recording had no valid pixels; likely too much motion')
    P = [];
    summaryEroded = nan(sz(1:2));
    return
end
if nargin<4
    doPlot = false;
end

if isempty(vIM) %for raster imaging
    vIM = ones(sz);
end

%initialize filtered image
IMf = IM; clear IM;
IMf(repmat(~valid, 1, 1, nTimePoints)) = nan;
nans = isnan(IMf);
nanFrac = mean(nans,3);
vIM(nans) = nan; %1000*mean(vIM(:,:, 1:min(end,400)), 'all', 'omitnan');

% %fill in missing values
% IMf= reshape(IMf, sz(1)*sz(2), []);
% incomplete = nanFrac>0 & nanFrac<1;
% IMs = nan(size(IMf));
% IMs(incomplete(:),:) = smoothdata(IMf(incomplete(:),:), 2, 'movmean', baselineWindow, 'omitnan');
% IMf(reshape(nans, size(IMf))) = IMs(reshape(nans, size(IMf))); clear IMs
% IMf = reshape(IMf, sz(1),sz(2), []);

if params.microscope == "SLAP2" || params.poissBasedStdIM
    %smooth the data at a timescale on which fluctuations look more
    %gaussian, for computing variances
    IMs = smoothdata(IMf./vIM, 3, 'movmean', ceil(denoiseWindow/2), 'omitnan');
    vIM = smoothdata(vIM, 3, 'movmean', ceil(denoiseWindow/2), 'omitnan');
    IMs = IMs.*vIM;

    %baseline estimate
    IMb = smoothdata(IMs, 3, 'movmedian', baselineWindow, 'omitnan');

    %estimate Vb and Vk, parameters for estimating variance from baseline brightness
    % Vb: the variance of a 'dim' pixel due to electronic and dark noise
    % Vk: the slope of the variance-brightness relationship
    firstValidFrames = find(any(~nans, [1 2]),500, 'first');
    varIM = var(IMs(:,:,firstValidFrames),0,3,"omitmissing");
    varIM(nanFrac>0.4) = nan;
    VIF=2;
    Vb = VIF*prctile(varIM, 10, 'all');
    varPred = mean(IMb(:,:,firstValidFrames),3,'omitmissing').* mean(vIM(:,:,firstValidFrames),3,'omitmissing');
    selBright = varPred>prctile(varPred(:), 90);
    Vk = prctile((varIM(selBright)-(Vb/VIF))./varPred(selBright), 10);

    %Highpass filter in time; This must occur before DoG to avoid edge artifacts
    IMf = IMf - IMb; 

    stdIM = sqrt(Vk.*IMb.*vIM+Vb); %compute standard deviation
else
    IMfden = smoothdata(IMf, 3, 'movmean', ceil(denoiseWindow/2), 'omitnan');
    %Highpass filter in time; This must occur before DoG to avoid edge artifacts
    IMb = smoothdata(IMfden, 3, 'movmedian', baselineWindow, 'omitnan');
    IMf = IMf - IMb;   %- smoothdata(IMf, 3, 'movmedian', baselineWindow, 'omitnan');

    % MAD-based robust standard deviation estimate
    stdIM = movmad(IMfden - IMb,baselineWindow,3,'omitmissing') ./ 0.6741891400433162.*ceil(denoiseWindow/2);
end
%divide by uncertainty to get a Z-score
IMf = IMf./stdIM;
clear IMb vIM

%time matched filter
gamma = exp(-1/tau);
mem = max(0,gamma*IMf(:,:,end));
for t = size(IMf,3):-1:1
    IMt = IMf(:,:,t);
    nanst = isnan(IMt);
    IMt(nanst) = mem(nanst);
    IMf(:,:,t) = gamma*mem + (1-gamma)*IMt;
    mem = IMf(:,:,t);
end
IMf(nans) = nan;

%Difference of Gaussians
IMf(nans) = 0;
IMf = imgaussfilt(IMf, [sigma sigma]);
IMf = IMf - imgaussfilt(IMf, 5*[sigma sigma]);
IMf(nans) = nan;

%nonmax suppression- find maxima
skIm = zeros(sz(1:2));
for fr = size(IMf,3)-ceil(1.5*tau):-1:2 %ceil(tau) because the filtering is uncertain in the final frames
    
    IMfr = IMf(:,:,fr); IMpre = IMf(:,:,fr-1); IMpost = IMf(:,:,fr+1);
    
    selMax = IMfr==ordfilt2(IMfr,9, ones(3));
    IMlocalMax(:,:,fr) = selMax & IMfr>IMpre & IMfr>=IMpost;
    %maxinds = find(IMfr==ordfilt2(IMfr,9, ones(3)));
    % sel = IMfr(maxinds)>0 & IMpre(maxinds)<=IMfr(maxinds) & IMpost(maxinds)<=IMfr(maxinds);
    % %sel = IMpre(maxinds)<=IMfr(maxinds) & IMpost(maxinds)<=IMfr(maxinds);
    % maxinds = maxinds(sel);
    skIm(IMlocalMax(:,:,fr)) = skIm(IMlocalMax(:,:,fr)) + IMfr(IMlocalMax(:,:,fr)).^2; 
end
skIm = skIm./(300+sum(~nans(:,:,2:end-ceil(1.5*tau)),3)); %normalize to # observations w regularizer
clear nans

%summary = skewness(IMf(:,:, 1:end-3*ceil(tau)), 1,3); %.*IMgamma; 
summaryEroded = skIm;

summaryEroded = summaryEroded ./ 10^(floor(log10(max(summaryEroded(:))))-1);
summaryEroded(~valid) = nan;
mfSummary = nanmedfilt2(summaryEroded, [5 5]);
summaryEroded = summaryEroded - mfSummary;
%valid = valid & (skIm ~= 0);
summaryEroded(~valid) = nan;

skIm(~valid) = nan;

thetaf = getActImPeaks(summaryEroded,params.peakth,[],params.peakFuncOpt,params.actImHeteroscedasticNoise,params.peakBufferSize);

P.row = thetaf(:,2);
P.col = thetaf(:,3);
P.val = thetaf(:,1);

peaksMask = zeros(size(summaryEroded));
peaksMask(sub2ind(size(peaksMask), round(P.row), round(P.col))) = 1;

P.peakIM = summaryEroded .* peaksMask;

if doPlot
    figure, imagesc(summaryEroded); %hAx1 = gca;
    hold on, scatter(P.col, P.row, 20*P.val, 'r' ); %'margeredgecolor', 'r');
    figure, imagesc(summaryEroded); %hAx2 = gca;
end

end

function val = gaussianPeaksIntegrated(theta, yxdata)
% Same interface as original:
% theta  : N×4  [A, mux, muy, sigma]
% yxdata : M×2  pixel centers [x, y]
% val    : M×1  integrated Gaussian over each pixel

x = yxdata(:,1);
y = yxdata(:,2);

% --- infer pixel size from grid ---
xu = unique(x);
yu = unique(y);

dx = median(diff(xu));
dy = median(diff(yu));

% build edges
xEdges = [xu(1)-dx/2; (xu(1:end-1)+xu(2:end))/2; xu(end)+dx/2];
yEdges = [yu(1)-dy/2; (yu(1:end-1)+yu(2:end))/2; yu(end)+dy/2];

W = numel(xu);
H = numel(yu);

% reshape index map
[~, xIdx] = ismember(x, xu);
[~, yIdx] = ismember(y, yu);

% --- parameters ---
A  = theta(:,1).';   % 1×N
mx = theta(:,2).';
my = theta(:,3).';
s  = max(theta(:,4).', eps);

c   = sqrt(pi/2);
rt2 = sqrt(2);

% --- integrate in x (W×N) ---
xL = xEdges(1:end-1);
xR = xEdges(2:end);
Ix = c .* s .* ( ...
    erf((xR - mx)./(rt2*s)) - ...
    erf((xL - mx)./(rt2*s)) );

% --- integrate in y (H×N) ---
yB = yEdges(1:end-1);
yT = yEdges(2:end);
Iy = c .* s .* ( ...
    erf((yT - my)./(rt2*s)) - ...
    erf((yB - my)./(rt2*s)) );

% --- combine to full image ---
img = (Iy .* A) * Ix.';   % H×W

% --- return in original point ordering ---
val = img(sub2ind([H W], yIdx, xIdx));
end


function val = gaussianPeaks(theta, yxdata)

y = yxdata(:,1);
x = yxdata(:,2);

A  = theta(:,1).';
my = theta(:,2).';
mx = theta(:,3).';
s  = theta(:,4).';

inv2s2 = 1 ./ (2 * s.^2);

dx = x - mx;
dy = y - my;

E = exp( -(dx.^2 + dy.^2) .* inv2s2 );
val = E * A.';
end

function loss = objFun(theta, yxdata, zdata, lambda)
loss = mean((gaussianPeaksIntegrated(theta,yxdata) - zdata).^2) + lambda * mean(abs(theta(:,1)));
end