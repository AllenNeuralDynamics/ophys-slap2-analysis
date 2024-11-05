function [F0] = computeF0(Fin, denoiseWindow, cutWindow, algoNum)
%Fin:   fluorescence traces, #timepoints is first dimension
%tau:   the timescale, in samples, of transients; algorithms are tested for
%       sparse spikes with time constant of tau samples
%cutWindow:     sustained transients longer than cutwindow tend to be absorbed
%               into F0
if nargin<4
    algoNum = 1;
end

algos = {@algo1, @algo2, @algo3};
algo = algos{algoNum};
warning('off');
F0 = algo(Fin, denoiseWindow, cutWindow);
warning('on');
end

function [F0] = algo1(Fin, denoiseWindow, hullWindow)
%algorithm 1:
%   use a median filter to remove noise and spikes
%   apply a rolling convex hull on the median filtered data
%   discard regions near nans that are more sensitive to noise, and replace with
%   extrapolated smooth F0
T = size(Fin,1);
hullWindow = min(hullWindow, floor(T/4));
deltaDes = max(4, (denoiseWindow/6)); % a safe spacing for computing F0, to save compute time
sampleTimes = round(linspace(1, T, ceil(T/deltaDes)+1));
nSampsInHull = ceil(hullWindow/deltaDes);

%median filter to reduce noise
F0 = medfilt1(Fin,denoiseWindow, [],1,'omitnan', 'truncate');

%minimum convex hull-like operation to remove effect of transients
origsz = size(F0);
F0 = reshape(F0, size(F0,1), []);
for cix = 1:size(F0,2)
    if all(isnan(F0(:,cix)))
        continue
    end
      
    for dix= nSampsInHull:-1:1
        xi =sampleTimes(dix:nSampsInHull:end);
        F00(:,dix) = interp1(xi, F0(xi,cix), sampleTimes);
    end
    FF = min(F00,[],2, 'omitnan'); 

    %convert to a smooth curve, discarding untrustworthy samples
    doubt = sum(~isnan(F00),2)<ceil(nSampsInHull/2);

        F2 = FF; 
        if sum(~doubt)>2
                F2(doubt) = nan;
        end
        fill = smoothdata(F2,1,'movmean', 2*ceil(nSampsInHull/2)+1, 'omitnan');
        F2(isnan(F2)) = fill(isnan(F2));
        F2= smoothdata(F2,1,'movmean', 2*ceil(nSampsInHull/2)+1, 'omitnan');
        F0(:,cix) = interp1(sampleTimes,F2,1:T,"pchip");
end
F0 = reshape(F0, origsz);
%F0(isnan(Fin)) = nan;
end



function F0 = algo2(Fin, denoiseWindow, hullWindow)
%algorithm 2:
%very fast
%median filter to denoise, then windowed minimum
%assumes first dimension is time

origsz = size(Fin);    
F0 = reshape(Fin, size(Fin,1), []);
nans = isnan(F0);
F0 = smoothdata(F0, 1,'movmedian',denoiseWindow, 'omitnan');
F0 = smoothdata(-imdilate(-F0, ones(hullWindow,1)), 1, 'movmean', hullWindow, 'omitnan');
F0(nans) = nan;
F0 = reshape(F0, origsz);
end


function F0 = algo3(Fin, denoiseWindow, hullWindow)
%algorithm 3
% no median; iteratively reweighted mean
%should be less biased
origsz = size(Fin);    
cutoff = 1.5;
Fin = reshape(Fin, size(Fin,1), []);
Fin = smoothdata(Fin, 1,'movmedian',10*denoiseWindow, 'omitnan');

F0 = smoothdata(Fin, 1,'movmedian',10*hullWindow, 'omitnan');
tmp = Fin-F0;
noise = std(tmp,0,1, 'omitnan');

for iter = 1:5
    setZero = abs(tmp./noise)>cutoff;
    tmp(setZero) = nan;
    noise = std(tmp, 0,1,"omitnan");
    F0 = F0 + smoothdata(tmp, 1,'movmean',hullWindow, 'omitnan');
    if iter<5
        tmp = Fin-F0;
    end
end
F0 = reshape(F0, origsz);
end
