function [summary, P, params] = localizeFlashesBergamo(IM, aData, params, doPlot)
%inputs:
%IM:        3D recording, X x Y x Time
%aData:     alignment metadata
tau = params.tau_s./(params.frametime*params.dsFac); %time constant in frames
params.tau_frames = tau;
sigma = params.sigma_px; %space constant in pixels

if nargin<4
    doPlot = false;
end

denoiseWindow = params.denoiseWindow_samps;
baselineWindow = ceil(params.baselineWindow_Glu_s/(params.frametime*params.dsFac));

%nanFrames = squeeze(all(isnan(IM),[1 2]));
nans = isnan(IM);
IMavg = mean(IM,3, 'omitnan');
IMgamma = sqrt(max(0,IMavg));

BG = prctile(IMavg(~isnan(IMavg)), 10);

%subtract background, relevant for Bergamo only (not SLAP2)
IMf = IM-BG;
IMavg = IMavg-BG;

IMf(nans) = 0;

%create a Difference of Gaussians filter (mean 0) that highights spots
DoGfilt = zeros(2*ceil(10*sigma)+1);
DoGfilt(ceil(end/2), ceil(end/2)) = 1;
DoGfilt = imgaussfilt(DoGfilt, [sigma sigma])-imgaussfilt(DoGfilt, 5*[sigma sigma]);
DoGfilt = DoGfilt - mean(DoGfilt(:));

%Apply spatial filter
IMf = convn(IMf,DoGfilt, 'same');
IMstruct = convn(max(0,IMavg), DoGfilt, 'same'); % filtered structural image, used for decorrelating motion

%Highpass filter in time
IMf(nans) = nan;
IMf = IMf - smoothdata(IMf, 3, 'movmedian', baselineWindow, 'omitnan'); 
nans = isnan(IMf);

%remove motion-associated variance
IMf = decorrelateMotion(IMf, IMstruct, aData, params);

% normalize by expected poisson noise
IMavg_nans = isnan(IMavg);
Pnoise = max(IMavg,0); %squared poisson noise
Mnoise = sqrt(...
    convn(Pnoise - circshift(Pnoise, [1 0]) , DoGfilt, 'same').^2 +...
    convn(Pnoise - circshift(Pnoise, [-1 0]) , DoGfilt, 'same').^2 +...
    convn(Pnoise - circshift(Pnoise, [0 1]) , DoGfilt, 'same').^2 +...
    convn(Pnoise - circshift(Pnoise, [0 -1]) , DoGfilt, 'same').^2);

IMnoise = sqrt(convn(Pnoise, abs(DoGfilt), 'same')); % total signal contributing to each measurement; sqrt(Mnoise) is inappropriate (should be no sqrt, with appropriate proportionality factor for photons in Pnoise) but keeps easy unit scaling
IMnoise = IMnoise + prctile(IMnoise(~IMavg_nans), 33) + 0.2*sqrt(Mnoise); % add a noise floor to the dimmer pixels
IMnoise(IMavg_nans) = nan;
IMf = IMf./IMnoise;

%temporal matched filter
IMf(nans) = 0;
mem = IMf(:,:,end);
gamma = exp(-1/tau);
for t = size(IMf,3):-1:1
    IMf(:,:,t) = max(0,gamma*mem) + (1-gamma)*IMf(:,:,t);
    mem = IMf(:,:,t);
end
IMf(nans) = nan;

%compute a summary image based on skewness
summary = mean(IMf.^3, 3, 'omitnan');
summary(imdilate(isnan(summary), ones(5))) = nan; %remove noisy edges

P = getTiledPeaks(IMf, IMavg, summary);

%plot a figure
if doPlot
    plotSummary(summary, IMgamma, P);
end
end

function plotSummary(summary, IMgamma, P)
red = summary./prctile(summary(~isnan(summary)), 99);
cyan = IMgamma./prctile(IMgamma(~isnan(IMgamma)),99);
figure, imshow(cat(3, red, cyan,cyan));
hold on,
scatter(P.col,P.row,100*(P.val./mean(P.val)).^2, 'm')
end

function P = getTiledPeaks(IM, IMavg, summary)
valid = ~imdilate(isnan(IM), ones(5)); 
peaks = valid;
for dim = 1:3
    peaks = peaks & IM > circshift(IM,1,dim) &  IM > circshift(IM,-1,dim);
end

linInds = find(peaks(:));
vvv = IM(linInds);
[rrr,ccc,ttt] = ind2sub(size(peaks), linInds);

%TO DO:
%We may want to discard localizations on dimmer pixels as inherently less
%likely, could use IMavg for this
%The selected region already selects for pixels with a large summary value
%vvv_ = vvv.*sqrt(IMavg(sub2ind(size(IMavg), rrr,ccc))); %signals normalized to image brightness

sz = size(IM);
tilesize = 64;
tilestartsR = 1:tilesize/2:(sz(1)-tilesize/4);
tileendsR = min(sz(1), tilestartsR+tilesize-1);
tilestartsC = 1:tilesize/2:(sz(2)-tilesize/4);
tileendsC = min(sz(2), tilestartsC+tilesize-1);

summaryVals = summary(sub2ind(size(summary), rrr, ccc));

keep = false(1,length(ttt)); %which events to keep
vNorm = zeros(length(ttt),1); %the event sizes, Z-scored

for rix = 1:length(tilestartsR)
    for cix = 1:length(tilestartsC)
        selStats = rrr>=tilestartsR(max(1,rix-1)) & rrr<=tileendsR(min(end,rix+1))  & ccc>=tilestartsC(max(1,cix-1)) & ccc<=tileendsC(min(end,cix+1));

        S = summary(tilestartsR(max(1,rix-1)):tileendsR(min(end,rix+1)), tilestartsC(max(1,cix-1)):tileendsC(min(end,cix+1)));
        Sp = prctile(S(:), [1 33]);
        Sthresh =  Sp(2) + 5*(Sp(2)-Sp(1));
        selS = summaryVals>Sthresh;

        vals = vvv(selStats & selS);
        ptile = prctile(vals, [1 50]);
        vals = 3*(vals-ptile(2))./(ptile(2) - ptile(1));
        thresh = ptile(2) + 4*(ptile(2) - ptile(1)); % threshold is 2*[98% confint], corresponding to an SNR of ~6
        selTile = rrr>=tilestartsR(rix) & rrr<=tileendsR(rix)  & ccc>=tilestartsC(cix) & ccc<=tileendsC(cix) & vvv>thresh & selS; %the events within this tile that should be kept
        keep(selTile) = true;
        vNorm(selStats & selS) = max(vNorm(selStats & selS), vals);
    end
end

rrr= rrr(keep);
ccc = ccc(keep);
ttt = ttt(keep);
vvv = vNorm(keep);

%upsample for superresolution
pC = []; pR = [];
for peakIx = length(ttt):-1:1
    R = IM(rrr(peakIx)+(-1:1), ccc(peakIx), ttt(peakIx));
    C = IM(rrr(peakIx), ccc(peakIx)+(-1:1), ttt(peakIx));

    ratioR = min(1e6,(R(2) - R(1))/(R(2) - R(3)));
    dR = (1-ratioR)/(1+ratioR)/2;
    pR(peakIx) = rrr(peakIx)-dR;

    ratioC = min(1e6,(C(2) - C(1))/(C(2) - C(3)));
    dC = (1-ratioC)/(1+ratioC)/2;
    pC(peakIx) = ccc(peakIx)-dC;
end

P.row = pR;
P.col = pC;
P.t = ttt;
P.val = vvv;

end

function IM = decorrelateMotion(IM, IMavg, aData,window)
nanFrames = squeeze(all(isnan(IM),[1 2]));
IMnan = isnan(IM);
sz = size(IM);
nPCs = 20;
goodPixels = mean(~IMnan,3)>0.9;
IM2 = IM;
IM2(IMnan) = 0;
IM2 = reshape(IM2(:,:,~nanFrames), sz(1)*sz(2), []);
IM2 = IM2- mean(IM2,2);
imageGrads = cat(3, IMavg-imtranslate(IMavg,[1 0]), IMavg-imtranslate(IMavg,[0 1]), IMavg-imtranslate(IMavg,[-1 0]), IMavg-imtranslate(IMavg,[0 -1]));
imageGrads = reshape(imageGrads, sz(1)*sz(2),[]);
[U,S,V] = svds(IM2(goodPixels,:),nPCs);
b = mvregress(imageGrads(goodPixels,:),U);
correction = (imageGrads*b)*S*V';
IM(:,:, ~nanFrames) = IM(:,:, ~nanFrames) - reshape(correction, sz(1), sz(2), []);
end

function P = motionPredictors (aData, window)
%predictors are
C = aData.motionDSc;
R = aData.motionDSr;
C = (C-mean(C))./sqrt(mean(C.^2,'all', "omitnan"));
R = (R-mean(R))./sqrt(mean(R.^2,"all", 'omitnan'));

dC = diff(C); dC = ([0 dC] + [dC 0])/2;
dR = diff(R); dR = ([0 dR] + [dR 0])/2;

P = cat(1,ones(size(C)), dC, dC.^2,dR, dR.^2);
%if nargin>1 && ~isempty(window)
Cs = C-smoothdata(C', window, 'movmedian')';
Rs = R-smoothdata(R', window, 'movmedian')';
P = cat(1, P, Cs, Cs.^2,Rs, Rs.^2, Cs.*Rs);
%end
if isfield(aData, 'motionDSz')
    keyboard
end
end
