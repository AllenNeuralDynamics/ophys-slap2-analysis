function summarize_LoCo(dr_or_pathToTrialTable, paramsIn)
%PARAMETER SETTING
if nargin>1
    if ischar(paramsIn)  % Parse JSON String to Structure
        paramsIn = jsondecode(paramsIn);
    end
    params = setParams('summarize_LoCo', paramsIn);
else
    params = setParams('summarize_LoCo');
end
if ~nargin
    [trialTablefn, dr] =  uigetfile('*.mat', 'Select a trialTable file', '*trialTable*.mat' );
else
    %parse dr
    %_or_pathToTrialTable
    if exist(dr_or_pathToTrialTable, 'dir')
        dr = dr_or_pathToTrialTable;
        trialTablefn = 'trialTable.mat';
    else
        [dr trialTablefn ext] = fileparts(dr_or_pathToTrialTable);
        trialTablefn = [trialTablefn ext];
    end
end

if params.makeJSON
    pythonenv_dir = uigetdir(getenv("USERPROFILE"),'Select Python Environment Directory');
    disp(['python env: ' pythonenv_dir])
end

params.startTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

copyReadDeleteScanImageTiff([]); %make sure we can use the function in parallel loops

%confirm that all files exist
[trialTable, keepTrials] = verifyFiles(trialTablefn, dr, params);
% for dmdIx = 1:numel(trialTable.refStack)
%     trialTable.refStack{dmdIx}.IM = []; %this uses a lot of memory and we won't need it
%end
nDMDs = size(trialTable.filename,1); %the trial table has size #DMDs x # trials; Bergamo is treated as '1 DMD'
nTrials = size(trialTable.filename,2);

%parameters that depend only on the microscope, hidden from GUI
switch params.microscope
    case 'SLAP2'
        trialTable.fnRaw = trialTable.filename;
    case 'bergamo'
        params.analyzeHz = nan;
end

disp(['## SUMMARIZING' newline 'Folder:'])
disp(dr)

savedr = [dr filesep 'ExperimentSummary'];
% on CodeOcean /data is read-only and we save to /results
is_CodeOcean = ~(getenv("CO_CPUS") == "");
if is_CodeOcean
    savedr = strrep(savedr, '/data', '/results');
end
if ~exist(savedr, 'dir')
    mkdir(savedr);
end
fnsave = [savedr filesep 'SummaryLoCo-' datestr(now, 'YYmmDD-HHMMSS') '.mat'];

%call up a GUI for the user to define Soma ROI and regions to exclude
if params.drawUserRois
    fnAnn = [dr filesep 'ANNOTATIONS.mat'];
    if exist(fnAnn, 'file')
        load(fnAnn, 'ROIs')
    else
        for DMDix = 1:nDMDs
            %load image data
            firstValidTrial = find(keepTrials(DMDix,:),1,"first");
            [~, fn, ext] = fileparts(trialTable.fnRegDS{DMDix,firstValidTrial});
            IM = copyReadDeleteScanImageTiff([dr filesep fn ext]);
            IM = squeeze(mean(IM,[3 4], 'omitnan'));
            hROIs(DMDix) = drawROIs(sqrt(max(0,IM)), dr, fn);
            ROIs(DMDix).dr = dr;
            ROIs(DMDix).fn = fn;
        end
        for DMDix = 1:nDMDs
            waitfor(hROIs(DMDix).hF);
            ROIs(DMDix).roiData = hROIs(DMDix).roiData;
        end
        save(fnAnn, 'ROIs'); clear hROIs;
    end
else
    ROIs = [];
end

%PROCESS DATA
for DMDix = nDMDs:-1:1
    %load some metadata
    firstValidTrial = find(keepTrials(DMDix,:),1,"first");
    fn = trialTable.fnAdata{DMDix,firstValidTrial};
    load([dr filesep fn], 'aData');
    numChannels = aData.numChannels;
    params.numChannels = numChannels;
    params.alignHz = aData.alignHz;
    if ~strcmpi(params.microscope, 'SLAP2')
        params.analyzeHz = 1/aData.frametime; %analyze conventional recordings at the acquisitoin framerate
    end
    if isfield(aData, 'Z')
        exptSummary.Z(DMDix) = aData.Z;
    else
        exptSummary.Z(DMDix) = nan;
        warning('Alignment data missing Z-plane, likely out of date!!')
    end
    clear aData

    %set up parallelization
    if params.nParallelWorkers>1
        p = gcp('nocreate');
        if isempty(p)
            poolsize = 0;
        else
            poolsize = p.NumWorkers;
        end
        dd = dir([dr filesep trialTable.fnRegDS{DMDix, firstValidTrial}]);
        try
            fileSize = dd.bytes;
        catch
            error(['Error loading registered tiff:' trialTable.fnRegDS{DMDix, firstValidTrial} '\n' 'Are paths in your trial table valid?']);
        end
        if ispc
            userMemInfo = memory;
            memAvailable = userMemInfo.MemAvailableAllArrays;
        else
            [~, result] = unix('grep MemAvailable /proc/meminfo | awk ''{print $2}''');
            memAvailable = str2double(result) * 1024;  % Convert KB to bytes
        end
        maxWorkers = min(size(trialTable.filename,2), floor(0.13*memAvailable/fileSize));
        nWorkers = min(params.nParallelWorkers, maxWorkers);
        
        if poolsize~=nWorkers ||  ~strcmpi(class(p), 'parallel.ProcessPool')
            delete(gcp('nocreate'));
            parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM 
        end
    else
        delete(gcp('nocreate'));
    end

    %Perform Localizations
    disp('Loading data and performing localizations...')
    mIM = cell(1, nTrials); aIM = cell(1,nTrials); alignData = cell(1, nTrials); peaks = cell(1, nTrials); discardFrames = cell(1,nTrials); %rawIMs = cell(1,nTrials)
    fns = trialTable.fnRegDS(DMDix, :);
    parfor trialIx = 1:nTrials
        if keepTrials(DMDix,trialIx)
            [~, mIM{trialIx}, aIM{trialIx}, alignData{trialIx}, peaks{trialIx}, discardFrames{trialIx}]= loadAndProcessTrialAsync(dr, fns{trialIx}, numChannels, params); %rawIMs{trialIx}
        end
    end
    %Assemble same-sized mean images from different-sized trial means
    szm1 = max(cellfun(@(x)size(x,1),mIM)); szm2 = max(cellfun(@(x)size(x,2), aIM));
    meanIM = nan(szm1,szm2,numChannels, nTrials); activIM = nan(szm1,szm2,1, nTrials);
    for trialIx = 1:nTrials
        tmp =  mIM{trialIx};
        meanIM(1:size(tmp,1),1:size(tmp,2),:,trialIx) = tmp;
        tmp =  aIM{trialIx};
        activIM(1:size(tmp,1),1:size(tmp,2),:,trialIx) = tmp;
    end
    params.sz = size(meanIM, [1 2]);

    %Make template
    disp('Making template for aligning across trials...')
    maxshift = 5;
    M = squeeze(sum(meanIM, 3));
    samples = find(keepTrials(DMDix,:)); samples = samples(unique(round(linspace(1,length(samples),20))));
    template = makeTemplateMultiRoi(M(:,:,samples), maxshift);

    %align all mean images to template
    disp('Aligning across trials...')
    meanAligned = [];
    actAligned = nan(size(meanIM,1), size(meanIM,2),1,nTrials);
    corrCoeff = nan(1,nTrials);
    motOutput = nan(2,nTrials);
    Mpad = nan([size(template) size(M,3)]);
    Mpad(maxshift+(1:size(M,1)), maxshift+(1:size(M,2)),:) = M;
    %clear M

    fillval = min(template(:),[], 'omitnan')-1;
    tFFT = fft2(max(template, fillval));
    for trialIx = nTrials:-1:1
        if ~keepTrials(DMDix,trialIx) || all(isnan(activIM(:,:,1,trialIx)), 'all')
            disp(['skipping trial, dmd:' int2str(trialIx) ' ' int2str(DMDix)])
            continue %skip
        end
        disp(['trial: ' int2str(trialIx)])
        
        output1 = dftregistration_clipped(tFFT, fft2(max(Mpad(:,:,trialIx), fillval)),1,80);
        mot1 = [-output1(3) -output1(4)]; %xcorr2_nans(Mpad(:,:,trialIx), template, [-output1(3) ; -output1(4)], maxshift);
        [motOutput(:,trialIx), corrCoeff(trialIx)] = xcorr2_nans(Mpad(:,:,trialIx), template, round(mot1'), maxshift);
        [rr,cc] = ndgrid(1:size(meanIM,1), 1:size(meanIM,2));

        for chIx = 1:size(meanIM,3)
            meanAligned(:,:,chIx,trialIx) = interp2(meanIM(:,:,chIx,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
        end
        actAligned(:,:,1,trialIx) = interp2(activIM(:,:,1,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
    end
    %clear Mpad activIM

    %identify outliers in alignment quality to determine valid trials
    ccf = corrCoeff;
    corrThresh = min(0.90, median(ccf, 'omitnan')-2*std(ccf, 'omitmissing'));
    actValidPix = squeeze(mean(~isnan(actAligned(:,:,1,:)), [1 2]));
    validTrials= find(ccf(:)>corrThresh & actValidPix(:)>mean(actValidPix)/2);
    exptSummary.meanIM{DMDix} = mean(meanAligned(:,:,:,validTrials),4, 'omitnan');
    actIM = prctile(actAligned(:,:,:,validTrials),80,4);  %mean(actAligned(:,:,:,validTrials), 4, 'omitnan');
    nanFrac = mean(isnan(actAligned(:,:,:,validTrials)), 4);
    actIM(nanFrac>0.6) = nan;
    exptSummary.actIM{DMDix} = actIM;

    %accumulate peaks, only from valid trials
    peaksCat = struct;
    for vTrialIx = 1:length(validTrials)
        trialIx = validTrials(vTrialIx);
        if ~isfield(peaksCat, 'row')
            peaksCat.row = peaks{trialIx}.row - motOutput(1,trialIx);
            peaksCat.col = peaks{trialIx}.col - motOutput(2,trialIx);
            peaksCat.val = peaks{trialIx}.val;
        else
            peaksCat.row = cat(1, peaksCat.row, peaks{trialIx}.row - motOutput(1,trialIx));
            peaksCat.col = cat(1, peaksCat.col, peaks{trialIx}.col - motOutput(2,trialIx));
            peaksCat.val = cat(1, peaksCat.val, peaks{trialIx}.val);
        end
    end

    %select sources
    %strategy 1: find peaks directly on aligned activity image
    actIM = mean(actAligned(:,:,:,validTrials), 4, 'includenan');
    actIM = actIM ./ 10^(floor(log10(max(actIM(:))))-1);
    nanFrac = mean(isnan(actAligned(:,:,:,validTrials)), 4);
    actIM(nanFrac>params.nanThresh) = nan;
    medIM = nanmedfilt2(actIM, 5.*[1 1]);
    actIM = actIM-medIM; %subtract a local baseline
    exptSummary.actIM{DMDix} = actIM;
    sz = size(actIM);

    %Mask out somata from activity image
    somaMask = false(size(actIM));
    if ~isempty(ROIs)
        for rix = 1:numel(ROIs(DMDix).roiData)
            if contains(upper(ROIs(DMDix).roiData{rix}.Label), 'SOMA')
                tmp = ROIs(DMDix).roiData{rix}.mask;
                somaMask(1:size(tmp,1), 1:size(tmp,2)) = somaMask(1:size(tmp,1), 1:size(tmp,2)) | tmp;
            end
        end
    end

    thetaf = getActImPeaks(actIM,params.peakth,medIM,2,somaMask);

    totalPix = sum(~isnan(actIM(:)) & ~somaMask(:));
    if totalPix == 0 | isempty(thetaf)
        k = 0;
    else
        sources.R = round(thetaf(:,2));
        sources.C = round(thetaf(:,3));
        sources.V = thetaf(:,1);
        k = length(sources.R);
    end

    %select regions near synapses, aligned across movies
    selPix = false([sz(1:2) k]);
    params.selRadius = ceil(2*params.dXY);
    for sourceIx = k:-1:1
        rr = round(sources.R(sourceIx));
        cc = round(sources.C(sourceIx));
        selPix(rr,cc,sourceIx) = true;
        selPix(:,:,sourceIx) = imdilate(selPix(:,:,sourceIx), strel('disk',params.selRadius));
    end
    pxAlwaysValid = mean(isnan(meanAligned(:,:,1,validTrials)),4)<params.nanThresh;
    selPix = selPix & repmat(pxAlwaysValid, 1, 1, k); %ADJUST SELECTED PIXELS NOT TO INCLUDE POORLY MEASURED PIXELS

    %prune any sources that got clipped by pixel selection process
    keepSources = sum(selPix, [1 2])>5;
    sources.R = sources.R(keepSources);
    sources.C = sources.C(keepSources);
    selPix = selPix(:,:,keepSources);
    disp(['Number of sources: ' int2str(sum(keepSources))]);
        
    %for each file, load high res data and refine
    params.tau_full=params.tau_s*params.analyzeHz;
    params = setParamsExtractTrial(params);
    
    if isempty(ROIs) || isempty(ROIs(DMDix))
        roiData =[];
    else
        roiData = ROIs(DMDix).roiData;
    end

    if any(keepSources)
        fns = trialTable.fnRaw(DMDix,:);
            if strcmpi(params.microscope, 'SLAP2')
                fls = trialTable.firstLine(DMDix,:);
                els = trialTable.lastLine(DMDix,:);
                E = processAllTrials_Async(dr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, motOutput, roiData, validTrials, params);
            else %BERGAMO
                fls = cell(1,numel(fns)); %first frame; leave empty for most uses
                els = cell(1,numel(fns)); %last frame; leave empty for most uses
                E = processAllTrials_Async(dr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, motOutput, roiData, validTrials, params);
            end

        %per-trial images
        exptSummary.E(:,DMDix) = E; %experiment data
    end
    exptSummary.selPix{DMDix} = any(selPix,3);
    exptSummary.aData(:,DMDix) = alignData;
    exptSummary.userROIs{DMDix} = roiData;
    exptSummary.peaks{DMDix}= peaks;
    exptSummary.perTrialMeanIMs{DMDix} = meanIM;
    exptSummary.perTrialMeanIMsAligned{DMDix} = meanAligned;
    exptSummary.perTrialActIms{DMDix} = actIM;
    exptSummary.perTrialActIMsAligned{DMDix} = actAligned;
    exptSummary.perTrialAlignmentOffsets{DMDix} = motOutput; %the alignment vector for each trial

    clear meanAligned meanIM actAligned F0selDS E
end

params.endTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

%prepare file for saving
exptSummary.params = params;
exptSummary.trialTable = trialTable;
exptSummary.dr = dr;

%save
save(fnsave, 'exptSummary', "-v7.3");

if params.makeJSON
    try
        setenv("PYTHONHOME",pythonenv_dir)
        pyrunfile([fullfile(fileparts(mfilename('fullpath')), 'generate_processing_json_SLAP2_multiROI_raster.py') ' --mat_path "' fnsave '" --output_dir "' savedr '"']);
        disp('Saved processing.json')
    catch
        disp('Did not save processing.json')
    end
else
    disp('Did not save processing.json')
end

disp('Done summarize_LoCo')
end
