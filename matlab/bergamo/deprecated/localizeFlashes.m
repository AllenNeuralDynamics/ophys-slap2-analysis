function [summary, P] = localizeFlashes(IMf)
tau = 2;
sigma = 1.5;
bfwin = 2*ceil(20*tau)+1;

nans = isnan(IMf);
nansC = convn(nans, true(1,1,7), 'same')>0;
IMf(nansC) = nan;
IMavg = mean(IMf,3, 'omitnan');
IMgamma = sqrt(IMavg);

BG = prctile(IMavg(~isnan(IMavg)), 10);

%smooth photons in space
IMf = IMf-BG;
IMf(nansC) = 0;
IMf = imgaussfilt(IMf, [sigma sigma])./(imgaussfilt(single(~nansC), [sigma sigma])+0.5); %weight by number of measurements; The samples at the edges will be noisier

%subtract baseline
IMf = IMf - smoothdata(IMf,3,'movmedian',bfwin, 'omitnan');

% normalize by expected poisson noise
IMavg_nans = isnan(IMavg);
IMg = IMgamma; IMg(IMavg_nans)=0;
IMg = imgaussfilt(IMg, [sigma sigma])./(imgaussfilt(single(~IMavg_nans), [sigma sigma])+0.5);
IMg(IMavg_nans) = nan;
IMf = IMf./(max(IMg, prctile(IMg(~IMavg_nans), 50)));

%temporal filter
mem = zeros(size(IMf, [1 2]));
gamma = exp(-1/tau);
for t = 1:size(IMf,3)
    IMf(:,:,t) = max(0,gamma*mem) + (1-gamma)*IMf(:,:,t);
    mem = IMf(:,:,t);
end
IMf(nansC) = nan;

%select regions with the most activity to search in
summary = mean(IMf.^3, 3, 'omitnan');
summary(imdilate(isnan(summary), ones(5))) = nan; %remove noisy edges

validArea = ~imdilate(nansC, ones(5));
P = getTiledPeaks(IMf, validArea, summary);

%plot a figure
% red = summary./prctile(summary(:), 99);
% cyan = IMgamma./prctile(IMgamma(:),99);
% figure, imshow(cat(3, red, cyan,cyan));
% hold on,
% scatter(P.col,P.row,100*(P.val./mean(P.val)).^2, 'm')

end

function P = getTiledPeaks(IM, valid, summary)
peaks = valid;
for dim = 1:3
    peaks = peaks & IM > circshift(IM,1,dim) &  IM > circshift(IM,-1,dim);
end

linInds = find(peaks(:));
vvv = IM(linInds);
[rrr,ccc,ttt] = ind2sub(size(peaks), linInds);

sz = size(IM);
tilesize = 64;
tilestartsR = 1:tilesize/2:(sz(1)-tilesize/4);
tileendsR = min(sz(1), tilestartsR+tilesize-1);
tilestartsC = 1:tilesize/2:(sz(2)-tilesize/4);
tileendsC = min(sz(2), tilestartsC+tilesize-1);

summaryVals = summary(sub2ind(size(summary), rrr, ccc));

keep = false(1,length(ttt));

for rix = 1:length(tilestartsR)
    for cix = 1:length(tilestartsC)
        selStats = rrr>=tilestartsR(max(1,rix-1)) & rrr<=tileendsR(min(end,rix+1))  & ccc>=tilestartsC(max(1,cix-1)) & ccc<=tileendsC(min(end,cix+1));
        
        S = summary(tilestartsR(max(1,rix-1)):tileendsR(min(end,rix+1)), tilestartsC(max(1,cix-1)):tileendsC(min(end,cix+1)));
        Sthresh = prctile(S(:), 50);
        selS = summaryVals>Sthresh;

        vals = vvv(selStats & selS);
        ptile = prctile(vals, [1 50]);
        thresh = ptile(2) + 2*(ptile(2) - ptile(1)); %threshold is an SNR of ~6
        th(rix,cix) = thresh;
        selTile = rrr>=tilestartsR(rix) & rrr<=tileendsR(rix)  & ccc>=tilestartsC(cix) & ccc<=tileendsC(cix) & vvv>thresh & selS;
        keep(selTile) = true;
    end
end

rrr= rrr(keep);
ccc = ccc(keep);
ttt = ttt(keep);
vvv = vvv(keep);

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
