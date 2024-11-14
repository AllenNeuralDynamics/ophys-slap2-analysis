function analyzeBergamo_Peaks

[fnsave, drsave] = uiputfile([obj.dr filesep obj.fn(1:end-4) '_DMD' int2str(obj.selDMD) 'analyzed.mat']);

maskData = obj.roiMasks;

%for each valid trial, load the fast data
validTrials = find(obj.exptSummary.validTrials);
sData.traces = nan(1,1,1,1,1); %time, roi, channel, traceType, trial
sData.noise = nan(1,1,1,1); %roi, trial, channel, traceType
sData.motionC = nan(1,1); %time, roi, channel, traceType, trial
sData.motionR = nan(1,1); %time, roi, channel, traceType, trial
sData.alignError = nan(1,1); %time, roi, channel, traceType, trial

for trialIx = 1:length(validTrials) %TMP
    fnTrial = obj.exptSummary.fns{validTrials(trialIx)};
    dmdNumStringIndex = strfind(fnTrial, 'DMD1')+3;
    fnTrial(dmdNumStringIndex) = int2str(obj.selDMD);

    ind =strfind(fnTrial, '_REGISTERED');
    datFn = [fnTrial(1:ind-1) '.dat'];

    sData.fns{trialIx} = datFn;

    if ~exist([obj.exptSummary.dr filesep datFn], 'file')
        msgbox(['Dat File: ' datFn ' not found. Please select location...']);
        [~, obj.exptSummary.dr] = uigetfile([obj.dr filesep '*.dat']); %update directory
    end

    if trialIx==1
        meta = loadMetadata([obj.exptSummary.dr datFn]);
        linerateHz = 1/meta.linePeriod_s;
        dt = ceil(linerateHz/obj.analyzeHz);
    end

    %load metadata
    S2data = slap2.Slap2DataFile([obj.exptSummary.dr datFn]);
    numLines = S2data.totalNumLines;
    numFrames = floor(numLines/dt);
    frameIxs = dt*(1:numFrames);

    %load alignment data
    aData = loadAlignmentData([obj.exptSummary.dr datFn]);
    motionC = interp1(aData.DSframes,aData.motionDSc, frameIxs, 'linear', 'extrap');
    motionR = interp1(aData.DSframes,aData.motionDSr, frameIxs, 'linear', 'extrap');
    aError =  interp1(aData.DSframes,aData.aError, frameIxs, 'linear', 'extrap');
    motPreds = [motionC ; motionR ; motionC.^2 ; motionR.^2 ; motionC.^3 ; motionR.^3 ; motionC.*motionR ; ones(size(motionC))];

    %load the high time resolution data
    %extize
    for rix =maskData.nROIs:-1:1
        D{rix} = nan([maskData.bbox{rix}(2,2)-maskData.bbox{rix}(1,2)+1, maskData.bbox{rix}(2,1)-maskData.bbox{rix}(1,1)+1, obj.numChannels, numFrames]);
    end

    disp('Reading High-res data')
    for fix = 1:numFrames
        if mod(fix, 100)==1
            disp(['Frame ' int2str(fix) ' of ' int2str(numFrames) '...'])
        end
        for cix = 1:obj.numChannels
            fdata = S2data.getImage(cix, ceil(fix*dt), 2*dt, 1); %Channel, Frame, Timewindow, Zindex
            for rix = 1:maskData.nROIs
                %query rows
                Rq = aData.cropRow+motionR(fix)-1+(maskData.bbox{rix}(1,2):maskData.bbox{rix}(2,2));
                %query columns
                Cq = aData.cropCol+motionC(fix)-1+(maskData.bbox{rix}(1,1):maskData.bbox{rix}(2,1));
                %meshgrid uses X,Y convention not R,C
                [Cq,Rq] = meshgrid(Cq,Rq);
                D{rix}(:,:,cix,fix) = interp2(fdata,Cq,Rq);
            end
        end
    end

    [traces, noise] = obj.extractTraces(D, maskData, motPreds); %

    %traces dimensions: [time, roi, channel, traceType]
    sData.trialLength(trialIx) = size(traces,1);
    if size(traces,1)>size(sData.traces,1) %The first dimension needs to be extended with nans
        sData.traces(end+1:size(traces,1),:,:,:,:) = nan;
        sData.motionC(end+1:size(traces,1),:) = nan;
        sData.motionR(end+1:size(traces,1),:) = nan;
        sData.alignError(end+1:size(traces,1),:) = nan;
    end
    sData.traces(:,:,:,:,trialIx) = nan; %prefill with nans
    sData.motionC(:,trialIx) = nan;
    sData.motionR(:,trialIx) = nan;
    sData.alignError(:,trialIx) = nan;

    %accumulate per-trace data
    sData.traces(1:size(traces,1),1:length(D),1:obj.numChannels, 1:size(traces,4), trialIx) = traces; %time, roi, channel, traceType, trial
    sData.noise(1:length(D), 1:obj.numChannels, 1:size(traces,4), trialIx) = noise; %roi, channel, traceType, trial
    sData.motionC(1:length(motionC),trialIx) =  motionC;
    sData.motionR(1:length(motionC),trialIx) =  motionR;
    sData.alignError(1:length(motionC),trialIx) = aError;

end
sData.maskData = maskData;
sData.frametime = 1/obj.analyzeHz;
save([drsave filesep fnsave], 'sData')

%save figure
disp('saving figure')
exportgraphics(obj.hAx, [drsave filesep fnsave(1:end-4) '_RoisFigure.pdf'], 'ContentType', 'vector');

%cluster the localizations
[density peaks sourceR sourceC] = clusterLocalizations(P, params);

extractSources(IM, peaks, sourceR, sourceC,  params)

end