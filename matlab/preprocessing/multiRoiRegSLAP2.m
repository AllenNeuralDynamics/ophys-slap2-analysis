function multiRoiRegSLAP2(fullPathToTrialTable, paramsIn)

if ~nargin
    [fn, dr] = uigetfile('*trialTable*.mat');
else
    [dr, fn, ext] = fileparts(fullPathToTrialTable); fn = [fn ext]; 
end

%PARAMETER SETTING
if nargin>1
    params = setParams('multiRoiRegSLAP2', paramsIn);
else
    params = setParams('multiRoiRegSLAP2');
end
    
%load the trial Table, which sets correspondences between the two DMDs
load([dr filesep fn], 'trialTable');

%set up parallelization
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
nWorkers = min(params.nWorkers, numel(trialTable.filename));
if poolsize<nWorkers
    delete(gcp('nocreate'));
    if params.nWorkers<15
        warning('You are using few parallel workers!');
    end
    parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM %4-30-24, lowering processes again to prevent another error (18 --> 15)
end
nDMDs = size(trialTable.filename,1);
[dixs,fixs] = ndgrid(1:nDMDs,1:length(trialTable.trueTrialIx));
fnRegDS = cell(nDMDs,length(trialTable.trueTrialIx));
fnAdata = cell(nDMDs,length(trialTable.trueTrialIx));
parfor p_ix = 1:numel(fixs)
    f_ix = fixs(p_ix); DMD_ix = dixs(p_ix);
    [fnRegDS{p_ix}, fnAdata{p_ix}]= alignAsync(dr, trialTable, params, f_ix, DMD_ix);
end
trialTable.fnRegDS = fnRegDS;
trialTable.fnAdata = fnAdata;
trialTable.alignParams = params;
save([dr filesep fn], "trialTable")

disp('done multiRoiRegistration.')
end

function [fnwrite, fnAdata] = alignAsync(dr, trialTable, params, f_ix, DMD_ix);

fn = trialTable.filename{DMD_ix,f_ix};
fnW = ['E' int2str(trialTable.epoch(f_ix)) 'T' int2str(f_ix) 'DMD' int2str(DMD_ix)];
firstLine = trialTable.firstLine(DMD_ix,f_ix);
lastLine = trialTable.lastLine(DMD_ix, f_ix);
aData = params;

disp(['Aligning: ' [dr filesep fn]])
    fnwrite = [dr filesep fnW '_REGISTERED_DOWNSAMPLED-' int2str(aData.alignHz) 'Hz.tif'];
    fnAdata = [dr filesep fnW '_ALIGNMENTDATA.mat'];

    if ~params.overwriteExisting && exist(fnAdata, 'file') && exist(fnwrite, 'file')
        disp([fn ' is already aligned; skipping' newline 'To force realign, pass TRUE as second argument']);
        return
    end

    retry=0;
    while retry<10
        try
            S2data = slap2.Slap2DataFile([dr filesep fn]);
            retry = 10;
        catch
            retry = retry+1
            continue
        end
    end
    meta = loadMetadata([dr filesep fn]);
    linerateHz = 1/meta.linePeriod_s;
    dt = linerateHz/aData.alignHz;
    numChannels = S2data.numChannels;
    %numLines = lastLine-firstLine+1;
    
    %sanity checks
    assert(length(S2data.hDataFile.fastZs)==1); %single plane acquisitions only

    %%%%%Make an initial template
    %crosscorrelate each initial frame to each other
    disp('generating template')
    initFrames = ceil(firstLine+dt : dt : min(lastLine-dt, firstLine+40*dt));
    nInitFrames = length(initFrames);

    if nInitFrames==0
        disp(['File ' fn ' was very short! Skipping alignment']);
        return
    end
    for fix = nInitFrames:-1:1
        for cix = numChannels:-1:1
            Y(:,:,cix,fix) = S2data.getImage(cix, initFrames(fix), ceil(dt), 1);
        end
    end

    Y = squeeze(sum(Y,3));
    
    %make data smaller for alignment 
    trimRows = find(~all(isnan(Y(:,:,1)),2), 1, 'first'):find(~all(isnan(Y(:,:,1)),2), 1, 'last'); 
    trimCols = find(~all(isnan(Y(:,:,1)),1), 1, 'first'):find(~all(isnan(Y(:,:,1)),1), 1, 'last'); 
    Y = Y(trimRows, trimCols,:);
    sz = size(Y);

    R = ones(nInitFrames);
    motion = zeros(2,nInitFrames,nInitFrames);
    for f1 = 1:nInitFrames
        for f2 = (f1+1):nInitFrames
            [motion(:,f1,f2), R(f1,f2)] = xcorr2_nans(Y(:,:,f2), Y(:,:,f1), [0 ; 0], 3);
            motion(:,f2,f1) = -motion(:,f1,f2);
            R(f2,f1) = R(f1,f2);
        end
    end
    [bestR, maxind] = max(median(R));
    frameInds = find(R(:,maxind)>=bestR);
    
    [viewR, viewC] = ndgrid((1:(sz(1)+2*aData.maxshift))-aData.maxshift, (1:(sz(2)+2*aData.maxshift))-aData.maxshift); %view matrices for interpolation
    template = nan(2*aData.maxshift+sz(1), 2*aData.maxshift+sz(2));
    for fix = 1:length(frameInds)
        template(:,:,fix) = interp2(1:sz(2), 1:sz(1), Y(:,:,frameInds(fix)),viewC-motion(2,frameInds(fix), maxind), viewR-motion(1,frameInds(fix), maxind), 'linear', nan);
    end
    template = mean(template, 3); 
    T0 = template;
    clear Y Yhp;

    initR = 0; initC = 0;
    DSframes = ceil(firstLine:dt:lastLine);
    nDSframes= length(DSframes); %number of downsampled frames
    
    %for spatial downsampling, used to calculate alignment quality and for quicker visualization
    dsTimes = 2;
    dsSz = floor(size(template)./(2^dsTimes));
    A_ds = nan([dsSz nDSframes]);

    motionDSr = nan(1,nDSframes); 
    motionDSc = nan(1,nDSframes); %matrices to store the inferred motion
    aErrorDS = nan(1,nDSframes); %alignment error output by dftregistration
    %aRankCorrDS = nan(1,nDSframes); %rank correlation, a better measure of alignment quality
    %recNegErr = nan(1,nDSframes); %rectified negative error, a better measure of alignment quality

    %output TIF
    pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM; 250nm
    fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);

    registrationFailed = false;
    disp('Registering:');
    try
    for DSframeIx = 1:nDSframes
        M1 = S2data.getImage(1, DSframes(DSframeIx), ceil(dt), 1); %moving image Ch1
        M1 = M1(trimRows, trimCols);
        if numChannels==2
            M2 =  S2data.getImage(2, DSframes(DSframeIx), ceil(dt), 1); %moving image Ch2
            M2 = M2(trimRows, trimCols);
            M = M1+M2;
        else
            M = M1;
        end
          
        if ~mod(DSframeIx, 1000)
            disp([int2str(DSframeIx) ' of ' int2str(nDSframes)]);
        end

        Ttmp = mean(cat(3, T0,template),3, 'omitnan');
        T = Ttmp(aData.maxshift-initR + (1:sz(1)), aData.maxshift-initC+(1:sz(2)));
        
        [motOutput, corrCoeff] = xcorr2_nans(M, T, [0 ; 0], aData.clipShift);
        motionDSr(DSframeIx) = initR+motOutput(1);
        motionDSc(DSframeIx) = initC+motOutput(2);
        aErrorDS(DSframeIx) = 1-corrCoeff^2;

        A1 = interp2(1:sz(2), 1:sz(1), M1,viewC+motionDSc(DSframeIx), viewR+motionDSr(DSframeIx), 'linear', nan);
        fTIF.WriteIMG(single(A1));
        if numChannels==2
            A2 = interp2(1:sz(2), 1:sz(1), M2,viewC+motionDSc(DSframeIx), viewR+motionDSr(DSframeIx), 'linear', nan);
            fTIF.WriteIMG(single(A2));
            A =A1+A2;
        else
            A = A1;
        end

        %downsample in space
        dsTmp = A;
        for dsIx = 1:dsTimes
            dsTmp = dsTmp(1:2:2*floor(end/2), 1:2:2*floor(end/2)) + dsTmp(1:2:2*floor(end/2), 2:2:2*floor(end/2)) + dsTmp(2:2:2*floor(end/2), 1:2:2*floor(end/2)) + dsTmp(2:2:2*floor(end/2), 2:2:2*floor(end/2));
        end
        A_ds(:,:,DSframeIx) = dsTmp;

        % Asmooth = imgaussfilt(A, [0.66 0.66]); %smooth to reduce variance effects of subframe interpolation
        % selCorr = sel & ~isnan(template);
        % aRankCorrDS(DSframeIx) = corr(Asmooth(selCorr), template(selCorr), 'type', 'Spearman');
        % recNegErr(DSframeIx) =  mean(min(0, Asmooth(selCorr)-template(selCorr)).^2)./mean(template(selCorr).^2);

        sel = ~isnan(A);
        nantmp = sel & isnan(template);
        template(nantmp) = A(nantmp);
        template(sel) = (1-aData.alpha)*template(sel) + aData.alpha*(A(sel));
        
        initR = round(motionDSr(DSframeIx));
        initC = round(motionDSc(DSframeIx));
    end
    catch ME
        disp(ME);
        registrationFailed = true;
    end

    fTIF.close;

    %Compute alignment error
    offset = 10;
    recNegErr = nan(1,nDSframes);
    nChunks = ceil(nDSframes./(aData.alignHz*10)); %align 10 seconds at a time
    chunkEdges = round(linspace(1, nDSframes+1, nChunks));
    for chunkIx = 1:length(chunkEdges)-1
        t_ixs = chunkEdges(chunkIx):chunkEdges(chunkIx+1)-1;
        template = median(A_ds(:,:,t_ixs), 3, 'omitmissing');
        nanFrac = mean(isnan(A_ds(:,:,t_ixs)), 3);
        template(nanFrac>0.2) = nan;
        template_gamma = sqrt(max(0,template)+offset);
        template = repmat(template, 1,1,length(t_ixs));
        template(isnan(A_ds(:,:,t_ixs))) = nan;
        recNegErr(1,t_ixs) = sqrt(squeeze(mean((max(0, (template-A_ds(:,:,t_ixs))./template_gamma).^2), [1 2], 'omitnan')./mean(max(0,(template./template_gamma).^2, 'includenan'), [1 2], 'omitnan')));
    end

    if std(motionDSc)>1.5 || std(motionDSr)>1.5
        registrationFailed = true;
    end
    if registrationFailed
        disp(['REGISTRATION ERROR OCCURRED FOR FILE: ' fn newline 'YOU MAY NEED TO QC THIS FILE!' newline 'CONTINUING...'])
        return
    end

    %save alignment metadata
    aData.numChannels = numChannels;
    aData.frametime = 1/aData.alignHz;
    aData.DSframes = DSframes;
    aData.motionDSc = motionDSc;
    aData.motionDSr = motionDSr;
    aData.aError = aErrorDS;
    %aData.aRankCorrDS = aRankCorrDS;
    aData.recNegErr = recNegErr;
    aData.cropRow = trimRows(1)-aData.maxshift; %offset to add to ROIs to index into original recording
    aData.cropCol = trimCols(1)-aData.maxshift; %offset to add to ROIs to index into original recording

    disp('Getting online motion correction offsets')
    [aData.onlineXshift, aData.onlineYshift, aData.onlineZshift] = getOnlineMotion(S2data.hDataFile, DSframes);

    %CONVERTING DATAFILE IMAGES INTO THE SAVED TIFF IMAGE SPACE:
    aData.trimRows = trimRows; %used to remap images from the datafile into the space of the saved tiffs
    aData.trimCols = trimCols;%used to remap images from the datafile into the space of the saved tiffs
    aData.viewC = viewC;%used to remap images from the datafile into the space of the saved tiffs
    aData.viewR = viewR;%used to remap images from the datafile into the space of the saved tiffs
    %EXAMPLE CODE
    %  Y = S2data.getImage(channel, lineTime, deltaTime, zPos);
    %  Ytrimmed = Y(aData.trimRows, aData.trimCols);
    %  sz = size(Ytrimmed);
    %  Yshifted = interp2(1:sz(2), 1:sz(1), Ytrimmed,aData.viewC+motionC, aData.viewR+motionR, 'linear', nan);

    save(fnAdata, 'aData');
end

function meta = loadMetadata(datFilename)
    ix = strfind(datFilename, 'DMD'+digitsPattern(1));
    metaFilename = [datFilename(1:ix+3) '.meta'];
    meta = load(metaFilename, '-mat');
end