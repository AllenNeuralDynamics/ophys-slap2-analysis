function [IMsel, IMselRaw, F0,  W, selNans] = prepareNMFproblem(IMsel, W0, F0selDS, params)
if ~params.nmfBackgroundComps
    F0 = nan(size(IMsel));
    szDS = size(F0selDS);
    xo = 1:szDS(2);
    xq = linspace(1, szDS(2), size(IMsel,2));
    for rix = 1:size(F0selDS, 1)
        F0(rix,:) = interp1(xo,F0selDS(rix,:), xq);
    end
    %baselineWindow = ceil(params.baselineWindow_Glu_s/params.frametime);
    %F0 = computeF0(IMsel', params.denoiseWindow_samps*params.dsFac+1, baselineWindow, 1)';
    IMsel = IMsel - F0;
    
    IMselRaw = IMsel;
    IMsel = matchedExpFilter(IMsel, params.tau_full);
    
    selNans = isnan(IMsel);
    IMsel(selNans) = 0;
    IMselRaw(selNans) = 0;

    %IMsel = IMsel- computeF0(IMsel', params.denoiseWindow_samps*params.dsFac+1,baselineWindow, 1)';
    W = W0;
else
    error('fitting global background components not implemented')
    % NOT IMPLEMENTED! 
    %GOALS
    % %compute a background (F0) per pixel
    % %NMF the background, and add those components to the NMF spatial components
    % %temporal filter the residual (dF)
    % %IMsel is the temporal filtered residual + background
    % F0sel = computeF0(IMsel', params.denoiseWindow*params.dsFac+1, params.baselineWindow*params.dsFac+1, 2)'; %algo2 doesn't underestimate F0 as badly in noisy conditions
    % dFsel = IMsel - F0sel;    
    % dFselTf = matchedExpFilter(dFsel, params.tau_full);
    % IMsel = dFselTf+F0sel;
    % 
    % selNans = isnan(IMsel);
    % nValid = sum(~selNans,2); %the number of values in each row
    % selRows = nValid>median(nValid)*0.9; %rows that we will use to estimate the background
    % 
    % IMBG = IMsel;
    % IMBG = IMBG - mean(IMBG,2, 'omitmissing');
    % IMBG(selNans) = 0;
    % [~,~,VV] = svds(IMBG(selRows,:), params.nmfBackgroundComps);
    % pred = cat(2, VV, ones(size(VV,1),1));
    % IMBG2 = smoothdata(IMBG, 2, 'movmedian', params.baselineWindow, 'omitmissing');
    % %for each row, fill nans with best fit
    % for rix = 1:size(IMBG,1)
    %     b = regress(IMBG2(rix, ~selNans(rix,:))', pred(~selNans(rix,:),:));
    %     fitVal = pred*b;
    %     IMBG(rix, selNans(rix,:)) = fitVal(selNans(rix,:));
    % end
    % 
    % %initialize components using PCA
    % tmp = IMsel - mean(IMsel,2, 'omitnan');
    % tmp(isnan(tmp)) = 0;
    % [U,~,VV] = svds(tmp,2);
    % clear tmp
    % UU = cat(2, U(:,1), -U(:,1), U(:,2), -U(:,2));
    % UU(UU<0) = 0;
    % UU = UU./sqrt(mean(UU.^2,1));
    % selU = mean(UU,1) > 0.66*mean(UU,'all');%select only components that are more spatially distributed
    % UU= UU(:,selU)+0.1;
    % 
    % 
    % IMBG =  smoothdata(IMBG,2,'movmean', ceil(params.baselineWindow*params.dsFac/2), 'omitnan');
    % opts = statset('MaxIter', 30,  'Display', 'final');
    % [Wbg,~] = nnmf(max(0,IMBG),size(UU,2), 'algorithm', 'mult', 'w0', UU, 'options', opts);
    % 
    % 
    % W = cat(2, W0, Wbg);
end
end