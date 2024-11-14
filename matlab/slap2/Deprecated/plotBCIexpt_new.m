function [pData, hF] = plotBCIexpt_new(exptSummary)

%outputs:
%hF are the figure handles
%pData is the data used in the plots

anParams.taskWindow_s = [0.5 20]; %window of time, in seconds, during which to average somatic output to get mean trial activity

%make plots
hF2 = plotAllTraces(exptSummary, anParams);
hF3 = plotTrialAverages(exptSummary, anParams);
hF1 = plotTotalInput(exptSummary, anParams);







end

function pData = mergeDMDs(sData)
for dmdNum = 1:length(sData)
    Tt = sData(dmdNum).traces;
    Nn = sData(dmdNum).noise;
    if dmdNum==1 %first DMD
        pData.F = Tt;
        pData.noise = Nn;
        pData.dmdIxs = ones(1,size(Tt,2));
        pData.names = sData(dmdNum).maskData.names;
    else
        if size(Tt,1)>size(pData.F,1)
            pData.F((size(pData.F,1)+1):size(Tt,1),:,:,:,:) = nan;
        else
            Tt((size(Tt,1)+1):size(pData.F,1),:,:,:,:) = nan;
        end
        pData.F= cat(2,pData.F,Tt);
        pData.noise = cat(1, pData.noise, Nn);
        pData.dmdIxs = cat(2, pData.dmdIxs, dmdNum*ones(1,size(Tt,2)));
        pData.names = cat(2, pData.names,  sData(dmdNum).maskData.names);
    end
end

fields2D = {'alignError', 'motionC', 'motionR'};
for fIx = 1:length(fields2D)
    pData.(fields2D{fIx}) = accumulate2Dfield(sData, fields2D{fIx});
end

%select frames to censor due to poor alignment
pData.motionCensor = censorMotion(pData);
pData.F(pData.motionCensor) = nan;

pData.denoiseWindow = [0.5/sData(1).frametime 0.5/sData(1).frametime];
pData.hullWindow = [0.8/sData(1).frametime 0.8/sData(1).frametime];

disp('Computing F0...')
pData.F0 = nan(size(pData.F));
for Ch = 1:2
    pData.F0(:,:,Ch,:,:) = computeF0(pData.F(:,:,Ch,:,:), pData.denoiseWindow(Ch), pData.hullWindow(Ch)); %allow F0 fitting to differ for different channels
end
pData.dF = pData.F-pData.F0;
pData.dFF = (pData.F-pData.F0)./max(0,pData.F0);

pData.tau(1) = 0.02/sData(1).frametime;%Glutamate time constant in samples
pData.tau(2) = 0.05/sData(1).frametime; %Calcium time constant in samples
end

function A = accumulate2Dfield(sData, fieldname)
for dmdNum = 1:length(sData)
    Tt = sData(dmdNum).(fieldname);
    if dmdNum==1 %first DMD
        A = Tt;
    else
        if size(Tt,1)>size(A,1)
            A((size(A,1)+1):size(Tt,1),:,:,:,:) = nan;
        else
            Tt((size(Tt,1)+1):size(A,1),:,:,:,:) = nan;
        end
        A= cat(3,A,Tt);
    end
end
end

function censoredFrames = censorMotion(pData)
%frames with sharp motion derivatives
dM= sqrt(mean(cat(4, abs(pData.motionC([2:end end],:,:) -  pData.motionC), abs(pData.motionC([1 1:end-1],:,:) - pData.motionC), ...
    abs(pData.motionR([2:end end],:,:) -  pData.motionR), abs(pData.motionR([1 1:end-1],:,:) - pData.motionR)).^2,4, 'omitnan'));
thresh = 0.05; %in pixels RMS

%frames with bad alignment errors
m = median(pData.alignError,1, 'omitnan');
e = max(0.02, 2*std(pData.alignError,0,1, 'omitnan')); % error threshold, in normalized error of dftRegistration

censoredFrames = dM>thresh & pData.alignError>(m+e);

%always censor frames in the wrong Z-plane
if isfield(pData, 'motionZ')
    censoredFrames = censoredFrames | pData.motionZ>1; %in microns
end
if size(censoredFrames,3)>1
    censoredFrames = cat(2, repmat(permute(censoredFrames(:,:,1), [1 3 4 5 2]), [1 sum(pData.dmdIxs==1)]), repmat(permute(censoredFrames(:,:,2), [1 3 4 5 2]), [1 sum(pData.dmdIxs==2)]));
else
    censoredFrames = repmat(permute(censoredFrames(:,:,1), [1 3 4 5 2]), [1 size(pData.F,2)]);
end
end

function hF = plotAllTraces(exptSummary, anParams)
outputCh = 2;
%nEpochs = max(exptSummary.trialTable.epoch);
colors = [0 1 0 ; 1 0 1];

dFFsom = accumSoma(exptSummary, outputCh); %channel

for spineCh = [1,2]
    [dF, Dixs] = accumTraces(exptSummary, 'dFraw', spineCh); %fieldname, channel

    %Figure with 2 subplots: 1 for soma, 1 for all synapses
        %selEpoch = exptSummary.trialTable.epoch==epoch;


        %denoised = accumTraces(exptSummary, 'denoised', spineCh); %fieldname, channel
        %dFnoise = estimatenoise(dF(:,1:3:end,:),2);

        D = dF; %shiftdim(sqrt(pData(epoch).noise(:,:,pData(epoch).traceType,:)), -1); %normalized to noise estimate

        %get trial lengths
        keepFrames = true(1, size(D,2)*size(D,3));
        for trialIx = size(D,3):-1:1
            tl = find(~all(isnan(D(:,:,trialIx)),1), 1, 'last');
            if isempty(tl)
                tl = 0;
            end
            trialLength(trialIx) = tl;
            keepFrames(((trialIx-1)*size(D,2)+trialLength(trialIx)+1):(trialIx*size(D,2))) = false;
        end
        trialbounds = cumsum(trialLength);

        hF(spineCh) = figure('name', ['All trials, Ch' int2str(spineCh)]);

        hAxSoma = subplot(10,1,1);
        tmp = reshape(dFFsom, 1, []);
        tmp(:, ~keepFrames) = []; %remove ends of trials
        T = (1:size(tmp,2))/exptSummary.params.analyzeHz;
        plot(T, tmp, 'color', colors(outputCh,:)); hold on;
        for tb = 1:length(trialbounds)
            plot((trialbounds(tb)/exptSummary.params.analyzeHz)*[1 1], ylim, 'k-');
        end
        set(hAxSoma, 'xtick', [])
        ylabel('Z-score')

        hAxSpines = subplot(10,1,2:9);
        tmp = reshape(D, size(D,1), []); %tmp is now 2D [ROIs allTime]
        tmp(:, ~keepFrames) = []; %remove ends of trials

        tmpSel = true(size(tmp));
        tmpErr = std(tmp,0,2, 'omitmissing');
        % for iter = 1:3 %normalize to std, ignoring outliers
        %     for rix = 1:size(tmp,1)
        %         tmpSel(rix,:) = tmp(rix,:)./tmpErr(rix)<4;
        %         tmpErr(rix) = std(tmp(rix, tmpSel(rix,:)));
        %     end
        % end
        tmp= tmp./tmpErr;

        T = (1:size(tmp,2))/exptSummary.params.analyzeHz;
        for traceIx = 1:size(tmp,1)
            offset = 10*traceIx;
            plot(T, tmp(traceIx,:)+offset, 'color', [0.8 0.8 0.8] ); hold on; %pData(epoch).colors(spineCh,:)
        end
        for traceIx = 1:size(tmp,1)
            offset = 10*traceIx;
            tmp2 = tmp(traceIx,:);
            tmp2(tmp2<(median(tmp2, 'omitmissing')+3)) = nan;
            plot(T, tmp2+offset, 'color', colors(spineCh,:), 'linewidth', 2);
        end
        for tb = 1:length(trialbounds)
            plot((trialbounds(tb)/exptSummary.params.analyzeHz)*[1 1], ylim, 'k-');
        end
        set(hAxSpines, 'ytick', [], 'xticklabel', [], 'tickdir', 'out', 'box', 'on', 'xgrid', 'on');
        ylabel('Z-score')
        xlabel('time (s)');
        set([hAxSoma hAxSpines], 'box', 'off')
        linkaxes([hAxSoma hAxSpines], 'x');

    %Plot motion
    %plot X,Y motion, Z motion? recNegErr?
    % tmp = cat(3, exptSummary.aData{:,}.motionC, pData(epoch).motionR);
    % tmp = reshape(tmp, [],size(tmp,3));
    % tmp(~keepFrames,:) = [];
    % tmp = tmp-mean(tmp,1, 'omitnan');
    % nans = any(isnan(tmp),2);
    % tmp(nans,:) = 0;
    % [UU,SS,~] = svds(tmp, 2); %PCA motion down to top 2 dimensions
    % motPCs = (UU.*diag(SS)')/sqrt(2);
    % motPCs(nans,:) = nan;
    %
    % hAxMotion = subplot(10,1,10);
    % plot(T, motPCs(:,1), 'color', 'b');
    % hold on,
    % plot(T, motPCs(:,2), 'color', 'r');
    % for tb = 1:length(trialbounds)
    %     plot(sData(1).frametime*trialbounds(tb)*[1 1], ylim, 'k-');
    % end
    %ylabel('Brain Movement')


end

% %export to json for collaborators
% jsonS.epoch = epoch;
% jsonS.channelNames = {'Glutamate', 'Calcium'};
% jsonS.dataDescription = 'delta Photons';
% jsonS.data = permute(D, [1 2 3 5 4]);
% jsonS.dataDimensions = {'time (frames at 100Hz)', 'ROI#', 'Channel', 'Trial'};
% jsonS.isSoma = isSoma;
% jsonS.motion = single(cat(3, pData(epoch).motionC, pData(epoch).motionR));
% jString = jsonencode(jsonS);
% fileID = fopen(['epoch' int2str(epoch) '.json'], 'w');
% fwrite(fileID, jString);
% fclose(fileID);
end

function [dFF, F, F0] = accumSoma(exptSummary, ch)

%get soma
for dix = 1:length(exptSummary.userROIs) % for each DMD
    for rix = 1:length(exptSummary.userROIs{dix})
        if strcmpi(exptSummary.userROIs{dix}{rix}.Label, 'soma')
            somaIxR = rix;
            somaIxD = dix;
        end
    end
end

F = [];
nTrials = size(exptSummary.E,1);
for tix = 1:nTrials
    if isempty(exptSummary.E{tix,somaIxD})
        F(1,:,tix) = nan;
    else
        tmp = exptSummary.E{tix,somaIxD}.ROIs.F(somaIxR,:,ch);
        if size(tmp,2)>size(F,2)
            F(1,(size(F,2)+1):size(tmp,2),:) = nan; %pad F
        elseif size(tmp,2)<size(F,2)
            tmp(1,(size(tmp,2)+1):size(F,2)) = nan; %pad tmp
        end
        F(1,:,tix) = tmp;
    end
end

F0 = permute(computeF0(permute(F, [2,3,1]), 32,301, 1), [3 1 2]);
dFF = (F-F0)./F0;
end

function [F, Dixs] = accumTraces(exptSummary, field, ch)
%outputs an array of size
% [#ROIs, #timepoiunts, #trials] for the given field name

nTrials = size(exptSummary.E,1);
%traces are aligned to trial start and padded with nans to give them a
%common size

%figure out how many ROIs are in each DMD
for dix = 1:2
    tix = find(~cellfun(@isempty,exptSummary.E(:,dix)),1, 'first');
    nSites(dix) = size(exptSummary.E{tix, dix}.F0,1); %#ok<AGROW>
end
Dixs = [ones(1,size(exptSummary.E{tix, 1}.F0,1)) 2*ones(1,size(exptSummary.E{tix, 2}.F0,1))];
F = [];

for tix = 1:nTrials
    %concatenate the two DMDs
    if isempty(exptSummary.E{tix,1}) && isempty(exptSummary.E{tix,2})
        tmp = nan(sum(nSites),1);
    elseif isempty(exptSummary.E{tix,1}) && ~isempty(exptSummary.E{tix,2})
        tmp2 = exptSummary.E{tix,2}.(field)(:,:,ch);
        tmp = cat(1, nan(nSites(1), size(tmp2,2)), tmp2);
    elseif ~isempty(exptSummary.E{tix,1}) && isempty(exptSummary.E{tix,2})
        tmp1 = exptSummary.E{tix,1}.(field)(:,:,ch);
        tmp = cat(1, tmp1, nan(nSites(2), size(tmp1,2)));
    else
        tmp = cat(1, exptSummary.E{tix,1}.(field)(:,:,ch), exptSummary.E{tix,2}.(field)(:,:,ch));
    end

    if size(tmp,2)>size(F,2)
        F(1:size(tmp,1),(size(F,2)+1):size(tmp,2),:) = nan; %pad F
    elseif size(tmp,2)<size(F,2)
        tmp(1:size(tmp,1),(size(tmp,2)+1):size(F,2)) = nan; %pad tmp
    end
    F(:,:,tix) = tmp; %#ok<AGROW>
end
end

function hF = plotTrialAverages(exptSummary, anParams)
outputCh = 2;
nEpochs = max(exptSummary.trialTable.epoch);
colors = [0 1 0 ; 1 0 1];

dF = accumTraces(exptSummary, 'dFraw', 1); %fieldname, channel
dFFsom = accumSoma(exptSummary, outputCh); %channel

for epoch = 1:nEpochs
    selEpoch = exptSummary.trialTable.epoch==epoch;
    output = dFFsom(:,:,selEpoch);
    tt = (1:size(dF,2))./exptSummary.params.analyzeHz; %time vector

    %mean and standard error of output signal over trials
    mO = mean(output,3, 'omitnan');
    eO = std(output,0,3, 'omitnan')./sqrt(sum(~isnan(output),3)-1);
    eO(isnan(eO)) = 0;

    hF = figure('name', ['Epoch ' int2str(epoch) ' Soma trial average']);
    hAxAvgOut(epoch) = axes; %#ok<AGROW>
    plot(tt, squeeze(output), 'color', colors(outputCh,:));
    hold on, plot(tt, mO, 'k', 'linewidth', 2);
    xlabel('time (s)');
    ylabel('dFF');
    title(['Somatic output all trials EPOCH' int2str(epoch)])

    selTime = tt>anParams.taskWindow_s(1) & tt<anParams.taskWindow_s(2); %time window for averaging dFF within trial
    normTime = tt<anParams.taskWindow_s(1);
    trialMeans{epoch} = squeeze(mean(output(1,selTime,:),2,'omitnan'));
    controlMeans{epoch} = squeeze(mean(output(1,normTime,:),2,'omitnan'));
end
linkaxes(hAxAvgOut);
set(hAxAvgOut, 'xlim', [0 20])


%did the somatic activity increase over epochs?
if nEpochs>1 %if multi-epoch
    alltrials = cell2mat(trialMeans');
    edges = cumsum(cellfun(@length, trialMeans));
    figure('name', ['Soma activity across trials']); plot(alltrials, 'o-');
    hold on,
    for E = 1:length(edges)
        plot(edges(E)*[1 1], ylim, 'color', [0.5 0.5 0.5], 'linewidth', 2)
    end
    xlabel('trial number');
    ylabel('mean activity');

    figure('name', ['Soma mean activity across epochs']),
    for epoch = 1:nEpochs
        scatter(epoch*ones(1,length(trialMeans{epoch})), trialMeans{epoch},'MarkerEdgeColor',[0.5 0.5 0.5]);
        hold on, plot(epoch+[-0.2 0.2], mean(trialMeans{epoch}, 'omitnan')*[1 1], 'k', 'linewidth', 2);

        if epoch<nEpochs
            %t-test
            [p,h] = ranksum(trialMeans{epoch}, trialMeans{epoch+1},'tail','left');
            disp(['Wilcoxon epoch' int2str(epoch+1) ':']);
            p
        end
    end
    set(gca, 'xlim', [0 nEpochs+1], 'xtick', 1:nEpochs);
    xlabel('Epoch');
    ylabel('mean activity')
end
% figure,
% shadedErrorPlot(tt,mO,eO,pData.colors(outputCh,:))



% N = shiftdim(sqrt(pData.noise(:,:,pData.traceType,:)), -1);
%
% tmpD = D(:,~isSoma, inputCh,1,:); %all input
% tmpN = N(:,~isSoma, inputCh,1,:);
% mO = mean(tmpD,5, 'omitnan');
% ste = std(tmpD,0,5, 'omitnan')./sqrt(sum(~isnan(tmpD),5)-1);
% ste(~isfinite(ste)) = 0;
% scale = mean(tmpN,5, 'omitnan');
% time = (0:size(tmpD,1)-1)*pData.frametime;
% hF = figure,
% shadedErrorPlot(time, mO./scale + (1:size(mO,2)), real(ste)./scale, pData.colors(inputCh,:))

end

function hF = plotPairwise(sData, pData)
tmp = Dnorm(:,~isSoma, inputCh,1,:);
tmp = reshape(permute(tmp, [2 1 5 4 3]), size(tmp,2), []); %tmp is now 2D [ROIs allTime]
tmp = tmp + 10*(0:(size(tmp,1)-1)); %offset for plotting
tmp(:, ~keepFrames) = []; %remove ends of trials
T = sData(1).frametime*(1:size(tmp,2));
plot(T, tmp, 'color', pData.colors(somaCh,:));
set(hAxSoma, 'xtick', [])
ylabel('Z-score')


hF = figure;

end

function hF = plotTotalInput(exptSummary, anParams)
%make a plot with the soma activity and the sum glutamate signal (in
%photons). Compute crosscorrelations between input and output
inputCh = 1;
outputCh = 2;
nEpochs = max(exptSummary.trialTable.epoch);
colors = [0 1 0 ; 1 0 1];

dF = accumTraces(exptSummary, 'dFraw', inputCh); %fieldname, channel
dFFsom = accumSoma(exptSummary, outputCh); %channel

for epoch = 1:nEpochs
    selEpoch = exptSummary.trialTable.epoch==epoch;
    nTrials = sum(selEpoch);

    input = dF(:,:,selEpoch);
    output = dFFsom(:,:,selEpoch);

    %get trial lengths
    keepFrames = true(1, size(input,2)*size(input,3));
    for trialIx = nTrials:-1:1
        tl = find(~all(isnan(input(:,:,trialIx)),1), 1, 'last');
        if isempty(tl)
            tl = 0;
        end
        trialLength(trialIx) = tl;
        keepFrames(((trialIx-1)*size(input,2)+trialLength(trialIx)+1):(trialIx*size(input,2))) = false;
    end

    %
    % tmpOut = D(:,isSoma,outputCh,:,:);
    % %tmpOut = tmpOut - mean(tmpOut, [1 5], 'omitnan');
    % output = sum(tmpOut, 2, 'omitnan');%in case there are multiple soma ROIs, sum them

    % %deconvolve output using OASIS
    % for trialIx = 1:size(output,5)
    %     df= output(:,:,:,:,trialIx);
    %     [c, s, options] = deconvolveCa(df, 'foopsi','type', 'ar1','pars', exp(-1/20), 'optimize_pars', false, 'optimize_b', false);
    %     %figure, plot(df); hold on, plot(c); hold on, plot(s)
    %     out_deconv(:,1,1,1,trialIx) = s;
    % end

    tmpIn= input;
    tmpIn = tmpIn - mean(tmpIn, 2, 'omitnan');
    sumInput = sum(tmpIn, 1, 'omitnan'); % omitting nans assumed mean dF=0 where not measured
    sumInput(all(isnan(tmpIn),1)) = nan;
    figure ('name', ['Sum Input, Epoch ' int2str(epoch)]),

    hAx(1) = subplot(2,1,1);
    plot(output(:), 'color', colors(outputCh,:));
    xlabel('Output; trials concatenated');
    hAx(2)= subplot(2,1,2);
    plot(sumInput(:), 'color', colors(inputCh,:));
    xlabel('Sum Input; trials concatenated');
    linkaxes(hAx, 'x')

    %%%
    [r_io, r_io_shuff, lags] = xcorr_vs_shuffle(squeeze(sumInput), squeeze(output), 5*exptSummary.params.analyzeHz);
    hF(1) = figure('Name',['Epoch ' int2str(epoch) ' Input-Output X-Corr']);
    plot(lags/exptSummary.params.analyzeHz, r_io, 'color', 'b', 'linewidth', 2);
    hold on, plot(lags/exptSummary.params.analyzeHz,  r_io_shuff, 'color', [0.5 0.5 0.5],'linewidth', 2);
    hold on, plot([0 0],  get(gca, 'ylim'), 'k:','linewidth', 2);
    xlabel('lag (s)');
    ylabel('Xcorrelation input vs output')
    legend({'', 'trial shuffle'})

    [r_ii, r_ii_shuff, lags] = xcorr_vs_shuffle(squeeze(sumInput), squeeze(sumInput), 5*exptSummary.params.analyzeHz);
    hF(2) = figure('name', ['Epoch ' int2str(epoch) ' Input-Input autoCorr']);
    plot(lags/exptSummary.params.analyzeHz, r_ii, 'color', colors(inputCh,:),'linewidth', 2);
    hold on, plot(lags/exptSummary.params.analyzeHz,  r_ii_shuff, 'color', [0.5 0.5 0.5],'linewidth', 2);
    hold on, plot([0 0],  get(gca, 'ylim'), 'k:','linewidth', 2);
    xlabel('lag (s)');
    ylabel('Input Autocorrelation')
    legend({'', 'trial shuffle'})

    [r_oo, r_oo_shuff, lags] = xcorr_vs_shuffle(squeeze(output), squeeze(output), 5*exptSummary.params.analyzeHz);
    hF(3) = figure('name', ['Epoch ' int2str(epoch) ' Output-Output autoCorr']);
    plot(lags/exptSummary.params.analyzeHz, r_oo, 'color',  colors(outputCh,:), 'linewidth', 2);
    hold on, plot(lags/exptSummary.params.analyzeHz,  r_oo_shuff,'color', [0.5 0.5 0.5], 'linewidth', 2);
    hold on, plot([0 0],  get(gca, 'ylim'), 'k:','linewidth', 2);
    xlabel('lag (s)');
    ylabel('Output Autocorrelation')
    legend({'', 'trial shuffle'})
end
end

function [r, rshuff, lags] = xcorr_vs_shuffle(d1,d2, window)
%computes cross-correlation of d1 and d2 across trials and compares to a
%null of shuffling with all non-synchronous trials
rii = nan(2*window+1,size(d1,2));
rij = nan(2*window+1,size(d1,2), size(d2,2));
d1(isnan(d1)) = 0;
d2(isnan(d2)) = 0;
for ix1  = 1:size(d1,2)
    [rii(:,ix1), lags] = xcorr(d1(:,ix1), d2(:,ix1), window);

    for ix2 = 1:size(d2,2)
        if ix2~=ix1
            rij(:,ix1,ix2) = xcorr(d1(:,ix1), d2(:,ix2), window);
        else
            rij(:,ix1,ix2) = nan;
        end
    end
end
r = mean(rii,2, 'omitnan');
rshuff= mean(rij, [2 3], 'omitnan');
end