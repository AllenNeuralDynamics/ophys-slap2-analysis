function integrationRegistration(fullPathToTrialTable, paramsIn)

if ~nargin
    [fn, dr] = uigetfile('*trialTable*.mat');
else
    [dr, fn, ext] = fileparts(fullPathToTrialTable); fn = [fn ext]; 
end

%PARAMETER SETTING
if nargin>1
    params = setParams('integrationRegistration', paramsIn);
else
    params = setParams('integrationRegistration');
end

params.startTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

%load the trial Table, which sets correspondences between the two DMDs
load([dr filesep 'trialTable.mat'], 'trialTable');

%% set up parallelization
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
nWorkers = min([params.nWorkers, numel(trialTable.filename), feature('numcores')]);
if poolsize~=nWorkers
    delete(gcp('nocreate'));
    if nWorkers<15
        warning('You are using few parallel workers!');
    end
    disp(['Parallel workers:' int2str(nWorkers)])
    parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM %4-30-24, lowering processes again to prevent another error (18 --> 15)
end

%% make look up table for each DMD

trialTable.lookupFile = [dr filesep 'integrationRegLookupTable.mat'];

% make lookup table if it doesn't exist
if ~exist(trialTable.lookupFile)
    lookupTable = makeRefLookupTable(dr,trialTable,params);

else % load user selected lookup table
    fprintf("Loading lookup table... "); tic;
    load(trialTable.lookupFile, 'lookupTable');
    fprintf('done, took %f sec\n',toc);
end

%% conduct alignment in parallel
nDMDs = size(trialTable.filename,1);

[dixs,fixs] = ndgrid(1:nDMDs,1:length(trialTable.trueTrialIx));
if params.saveTiffs; fnRegDS = cell(nDMDs,length(trialTable.trueTrialIx)); end
fnAdata = cell(nDMDs,length(trialTable.trueTrialIx));
parfor p_ix = 1:numel(fixs)
    f_ix = fixs(p_ix); DMD_ix = dixs(p_ix);
    [fnRegDS{p_ix}, fnAdata{p_ix}]= alignIntegrationAsync(dr, trialTable, lookupTable, params, f_ix, DMD_ix);
end

params.endTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

if params.saveTiffs; trialTable.fnRegDSInt = fnRegDS; end
trialTable.fnAdataInt = fnAdata;
trialTable.alignParamsInt = params;
save([dr filesep fn], "trialTable")

disp('done integrationRegistration.')
end

function lookupTable = makeRefLookupTable(dr, trialTable, params)
nDMDs = size(trialTable.filename,1);
for DMDix = nDMDs:-1:1
    [~,n] = fileparts(trialTable.filename{DMDix,1});
    n_base = regexprep(n,'-TRIAL[0-9]+','','ignorecase');
    n_base = regexprep(n_base,'-CYCLE[0-9]+','','ignorecase');
    metaDataFileName = fullfile(dr, [n_base '.meta']);
    mustBeFile(metaDataFileName);
    metaData = load(metaDataFileName, '-mat');

    list = dir([dr filesep '**' filesep '*DMD' int2str(DMDix) '_CONFIG2-REFERENCE*']);
    if isempty(list)
        list = dir([dr filesep '**' filesep '*DMD' int2str(DMDix) '*-REFERENCE*']);
    end
    ReferenceStack_ = slap2.gui.refstack.ReferenceStack.loadTif(fullfile(list(1).folder, list(1).name));
    numChannelsRefStack = numel(ReferenceStack_.data);
    for chIx = numChannelsRefStack:-1:1
        refStack(:,:,:,chIx) = permute(ReferenceStack_.data{chIx},[2 1 3]);
    end

    if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
        numLinesPerCycle = length(metaData.AcquisitionContainer.AcquisitionPlan.superPixelIDs);

        zs_ix = horzcat(metaData.AcquisitionContainer.AcquisitionPlan.activeZs{:});
        zs_ix = unique(zs_ix);
        zs = zeros(length(zs_ix),1);

        true_zs = cellfun(@(x) x.z, metaData.AcquisitionContainer.ROIs.rois);

        zPlanes = nan(1,numLinesPerCycle);
        for ix = 1:numLinesPerCycle
            if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan.activeZs{ix})
                zPlanes(ix) = metaData.AcquisitionContainer.AcquisitionPlan.activeZs{ix}(1);
            end
        end
        for ix = 1:length(zs_ix)
            if all(isnan(metaData.AcquisitionContainer.AcquisitionPlan.zTrajectory))
                zs(ix) = metaData.remoteFocusPosition_um;
            else
                zPlane_um = mean(metaData.AcquisitionContainer.AcquisitionPlan.zTrajectory(zPlanes == zs_ix(ix)));
                [~, zPlane_trueIx] = min(abs(true_zs - zPlane_um));
                zs(ix) = true_zs(zPlane_trueIx);
            end
        end
    else
        zs = metaData.AcquisitionContainer.ParsePlan.zs;
        numLinesPerCycle = length(metaData.AcquisitionContainer.ParsePlan.acqParsePlan);
    end

    dmdPixelsPerColumn = metaData.dmdPixelsPerColumn;
    dmdPixelsPerRow = metaData.dmdPixelsPerRow;
    numFastZs = length(zs);

    % get list of superpixels and extract data

    allSuperPixelIDs{DMDix} = [];

    fastZ2RefZ{DMDix} = zeros(numFastZs,1);
    middleRefZ = (max(ReferenceStack_.zs) + min(ReferenceStack_.zs))/2;
    middleZ = (max(zs) + min(zs))/2;
    for z = 1:numFastZs
        [~, ind] = min(abs((ReferenceStack_.zs - middleRefZ) - (zs(z) - middleZ)));
        fastZ2RefZ{DMDix}(z) = ind;
    end

    for lineSweepIdx = 1:numLinesPerCycle
        if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
            superPixIdxs = metaData.AcquisitionContainer.AcquisitionPlan.superPixelIDs{lineSweepIdx}';
        else
            superPixIdxs = metaData.AcquisitionContainer.ParsePlan.acqParsePlan(lineSweepIdx).superPixelID;
        end

        if numel(superPixIdxs) == 0; continue; end
        
        if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
            zIdx = metaData.AcquisitionContainer.AcquisitionPlan.activeZs{lineSweepIdx}(1);
        else
            zIdx = metaData.AcquisitionContainer.ParsePlan.acqParsePlan(lineSweepIdx).sliceIdx(1)+1;
        end

        % only take integration mode pixels
        if params.integrationOnly
            superPixIdxs(superPixIdxs <= dmdPixelsPerColumn*dmdPixelsPerRow) = [];
        end

        spIDs = superPixIdxs*100+zIdx; % add z plane to end of superpixel ID
        allSuperPixelIDs{DMDix} = [allSuperPixelIDs{DMDix}; spIDs]; % make list of unique superpixels across all Zs
    end

    [allSuperPixelIDs{DMDix}, ~, ic] = unique(allSuperPixelIDs{DMDix});
    spSampleCt = accumarray(ic, 1);

    fprintf("%d superpixels detected\n", length(allSuperPixelIDs{DMDix}));

    % make sparse matrix with each superpixel's corresponding mask (roiMasks)

    fprintf("Calculating ROI Masks... ");
    tic;

    % using sparse matrix
    sparseMaskInds{DMDix} = [];
    if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
        allPixelReplacementMaps = metaData.AcquisitionContainer.AcquisitionPlan.pixelReplacementMaps;
    else
        allPixelReplacementMaps = metaData.AcquisitionContainer.ParsePlan.pixelReplacementMaps;
    end

    for i = 1:length(allSuperPixelIDs{DMDix})
        tmpMask = zeros(dmdPixelsPerColumn,dmdPixelsPerRow,numFastZs);

        sp = allSuperPixelIDs{DMDix}(i);
        zIdx = rem(sp,100);
        superPixIdx = (sp - zIdx) / 100;
        pixelReplacementMap = allPixelReplacementMaps{zIdx};

        open = uint32(pixelReplacementMap(pixelReplacementMap(:,2) == superPixIdx,1))+1;

        if isempty(open)
            open = superPixIdx+1;
        end

        openR = idivide(open-1, dmdPixelsPerRow, 'floor')+1;
        openC = open - (openR-1) * dmdPixelsPerRow;
        openPixs = uint32(openR + (openC-1) * dmdPixelsPerColumn + double(zIdx-1) * dmdPixelsPerColumn * dmdPixelsPerRow);

        sparseMaskInds{DMDix} = [sparseMaskInds{DMDix}; openPixs, ones(size(openPixs))*i];
    end
    clear('tmpMask');

    roiMasks = sparse(sparseMaskInds{DMDix}(:,1),sparseMaskInds{DMDix}(:,2),1,dmdPixelsPerColumn * dmdPixelsPerRow * numFastZs,length(allSuperPixelIDs{DMDix}));
    fprintf('done DMD%d - took %f sec\n', DMDix, toc);

    % make sure refStack is larger than imaged frame

    bl = single(refStack)/100; %  reference image (X x Y x Z x C)
    bl(bl < 0) = 0; % remove any negative values, should not happen
    bl = padarray(bl,[floor((max(size(bl,1),dmdPixelsPerColumn)-size(bl,1))/2),...
        floor((max(size(bl,2),dmdPixelsPerRow)-size(bl,2))/2),...
        0],...
        mean(bl(:)),...
        'both');
    
    assert(size(bl,3) > numFastZs, 'ERROR: fewer Z planes in reference stack than fastZs');

    if size(bl,1) < dmdPixelsPerColumn
        bl = padarray(bl,[1,0,0],mean(bl(:)),'pre');
    end

    if size(bl,2) < dmdPixelsPerRow
        bl = padarray(bl,[0,1,0],mean(bl(:)),'pre');
    end

    % Calculate lookup table

    xPre = params.maxshiftXY; xPost = params.maxshiftXY;
    yPre = params.maxshiftXY; yPost = params.maxshiftXY;
    zPre{DMDix} = min(params.maxshiftZ, min(fastZ2RefZ{DMDix} - 1));
    zPost{DMDix} = min(params.maxshiftZ, min(size(bl,3) - fastZ2RefZ{DMDix}));

    likelihood_means{DMDix} = makeLookupTable(bl, sparseMaskInds{DMDix}, numFastZs, fastZ2RefZ{DMDix},-yPre:yPost,-xPre:xPost,-zPre{DMDix}:zPost{DMDix},ReferenceStack_.channels);
    % likelihood_means{DMDix} = likelihood_means{DMDix} .* repmat(reshape(spSampleCt,[1 1 1 1 length(allSuperPixelIDs{DMDix})]),[size(likelihood_means{DMDix},1:4) 1]);
end
lookupTable.likelihood_means = likelihood_means;
lookupTable.allSuperPixelIDs = allSuperPixelIDs;
lookupTable.sparseMaskInds = sparseMaskInds;
lookupTable.xPre = xPre;
lookupTable.xPost = xPost;
lookupTable.yPre = yPre;
lookupTable.yPost = yPost;
lookupTable.zPre = zPre;
lookupTable.zPost = zPost;
lookupTable.fastZ2RefZ = fastZ2RefZ;

save(trialTable.lookupFile,'lookupTable','-v7.3');

end


function [fnwrite, fnAdata] = alignIntegrationAsync(dr, trialTable, lookupTable, params, f_ix, DMD_ix)

fn = trialTable.filename{DMD_ix,f_ix};
fnW = ['E' int2str(trialTable.epoch(f_ix)) 'T' int2str(f_ix) 'DMD' int2str(DMD_ix) '_INTEGRATION'];
firstLine = trialTable.firstLine(DMD_ix,f_ix);
lastLine = trialTable.lastLine(DMD_ix, f_ix);
aData = params;

disp(['Aligning: ' [dr filesep fn]])

hSlap2DataFile = slap2.Slap2DataFile([dr filesep fn]);
if isprop(hSlap2DataFile, 'hDataFile')
    hLowLevelDataFile = hSlap2DataFile.hDataFile;
else
    hLowLevelDataFile = hSlap2DataFile.hMultiDataFiles;
end
numLinesPerCycle = hLowLevelDataFile.header.linesPerCycle;
totalCycles = hLowLevelDataFile.numCycles;
numChannels = hSlap2DataFile.numChannels;
channels = hLowLevelDataFile.metaData.channelsSave;

assert(length(channels) == numChannels, 'Saved channels does not match numChannels!');

linerateHz = 1/hLowLevelDataFile.metaData.linePeriod_s;
dt = linerateHz/aData.alignHz;

if dt < numLinesPerCycle
    aData.alignHz = floor(linerateHz / numLinesPerCycle);
    dt = linerateHz/aData.alignHz;
    disp(['Requested DS freq too fast, adjusting to ' int2str(aData.alignHz) ' Hz']);
end

if params.saveTiffs
    fnwrite = [dr filesep fnW '_REGISTERED_DOWNSAMPLED-' int2str(aData.alignHz) 'Hz.tif'];
    if params.efficientTiffSave
        fnwriteTmp = [params.tempFileDir filesep int2str(round(rand(1)*10000)) '.tif'];
    else
        fnwriteTmp = fnwrite;
    end
else
    fnwrite = 'no tiff saved';
end

fnAdata = [dr filesep fnW '_ALIGNMENTDATA.mat'];

if ~params.overwriteExisting && exist(fnAdata, 'file') % && exist(fnwrite, 'file')
    disp([fn ' is already aligned; skipping' newline 'To force realign, pass TRUE as second argument']);
    return
end

DSframes = ceil(firstLine:dt:lastLine);
nDSframes= length(DSframes); %number of downsampled frames

%% load metadata
dmdPixelsPerColumn = hLowLevelDataFile.metaData.dmdPixelsPerColumn;
dmdPixelsPerRow = hLowLevelDataFile.metaData.dmdPixelsPerRow;
numFastZs = length(hLowLevelDataFile.fastZs);

%%  Calculate metric table and infer motion from MAP
if isfield(aData, 'motionMetric')
    motionMetric = lower(string(aData.motionMetric));
    motionMetric = erase(motionMetric, "'");
else
    motionMetric = "poisson";
end
validMotionMetrics = ["poisson", "correlation"];
assert(any(motionMetric == validMotionMetrics), ['Unknown motionMetric: ' char(motionMetric)]);

% how much change is allowed in each dimension at each step
searchRadius = aData.clipShift;
searchRadiusZ = ceil(searchRadius / 2);

motionDS = nan(nDSframes,3);
brightnessDS = nan(nDSframes,1);
loglikelihoodDS = nan(nDSframes,1);

% dataMatrix = zeros(length(lookupTable.allSuperPixelIDs{DMD_ix}),nDSframes);
% expectedMatrix = zeros(length(lookupTable.allSuperPixelIDs{DMD_ix}),nDSframes);

if params.saveTiffs
    pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM; 250nm
    fTIF = Fast_BigTiff_Write(fnwriteTmp,pixelscale,0);
end

xMotRange = lookupTable.xPre + lookupTable.xPost + 1;
yMotRange = lookupTable.yPre + lookupTable.yPost + 1;
zMotRange = lookupTable.zPre{DMD_ix} + lookupTable.zPost{DMD_ix} + 1;

ySearch = 1:yMotRange;
xSearch = 1:xMotRange;
zSearch = 1:zMotRange;

% Get unique superpixel IDs and their group memberships
[~, ~, groupIndices] = unique(lookupTable.sparseMaskInds{DMD_ix}(:, 2));

% Group the pixel indices by superpixel
splitPixels = accumarray(groupIndices, lookupTable.sparseMaskInds{DMD_ix}(:, 1), [], @(x) {x});

% Compute the median pixel index for each superpixel
medianIndices = cellfun(@(x) x(round(length(x)/2)), splitPixels);

% Convert linear indices to row and column coordinates
[spRows, spCols, spPlanes] = ind2sub([dmdPixelsPerColumn, dmdPixelsPerRow, numFastZs], medianIndices);

registrationFailed = false;
fprintf("Inferring motion...\n")
tic
try
for DSframeIx = 1:nDSframes
    timeWindow = max(1,floor(DSframes(DSframeIx)-2*dt)):min(ceil(DSframes(DSframeIx)+2*dt),hLowLevelDataFile.totalNumLines);

    if ~mod(DSframeIx, 1000)
        disp([int2str(DSframeIx) ' of ' int2str(nDSframes)]);
    end

    lineIndices  = mod(timeWindow-1,numLinesPerCycle)+1;
    cycleIndices = floor((timeWindow-1) / numLinesPerCycle)+1;

    % load data
    data = zeros(length(lookupTable.allSuperPixelIDs{DMD_ix}),numChannels);
    spCt = zeros(length(lookupTable.allSuperPixelIDs{DMD_ix}),1);
    allLineData = hLowLevelDataFile.getLineData(lineIndices, cycleIndices, channels);
    for t = 1:length(allLineData)
        superPixIdxs = hLowLevelDataFile.lineSuperPixelIDs{lineIndices(t)};

        if numel(superPixIdxs) == 0; continue; end

        lineData = max(0,allLineData{t});
        zIdx = hLowLevelDataFile.lineFastZIdxs(lineIndices(t));

        spID = superPixIdxs*100 + uint32(zIdx); % make superpixel index with Z plane
        [~,spIdx] = ismember(spID,lookupTable.allSuperPixelIDs{DMD_ix});

        data(spIdx(spIdx>0),:) = data(spIdx(spIdx>0),:) + single(lineData(spIdx>0,:));
        spCt(spIdx(spIdx>0)) = spCt(spIdx(spIdx>0)) + 1;
    end
    data = data ./ 100; %./ spSampleCt;
    data(spCt == 0,:) = nan;

    if mean(spCt == 0) > 0.5; continue; end

    % calculate motion score at all shifts
    expectedMeans = bsxfun(@times, lookupTable.likelihood_means{DMD_ix}, reshape(spCt, [1 1 1 1 length(lookupTable.allSuperPixelIDs{DMD_ix})]));
    if motionMetric == "poisson"
        % log of mean before brightness scale s must use same μ_ref·spCt as expectedMeans (Poisson k·log λ with λ = s·μ_ref·spCt)
        log_means = log(expectedMeans + 1e-8);
        [logLikelihoodTable, scalingFactorTable] = poissonLogLikelihoodTable(data, expectedMeans, ...
                                                                            log_means,ySearch,xSearch,zSearch,channels,params.robust);
    else
        [logLikelihoodTable, scalingFactorTable] = correlationLikelihoodTable(data, expectedMeans, ySearch, xSearch, zSearch, channels);
    end

    [loglikelihoodDS(DSframeIx), I] = max(logLikelihoodTable(:));

    [My, Mx, Mz] = ind2sub(size(logLikelihoodTable),I);

    %perform superresolution upsampling, assuming quadratic loglikelihood
    dY = 0;dX = 0;dZ = 0;
    if My>1 && My<size(logLikelihoodTable,1)
        ratioY = min(1e6,(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My-1,Mx,Mz))/(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My+1,Mx,Mz)));
        dY = (1-ratioY)/(1+ratioY)/2;
    end
    if Mx>1 && Mx<size(logLikelihoodTable,2)
        ratioX =min(1e6, (logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx-1,Mz))/(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx+1,Mz)));
        dX = (1-ratioX)/(1+ratioX)/2;
    end
    if Mz>1 && Mz<size(logLikelihoodTable,3)
        ratioZ =min(1e6, (logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx,Mz-1))/(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx,Mz+1)));
        dZ = (1-ratioZ)/(1+ratioZ)/2;
    end

    motionDS(DSframeIx,:) = [ySearch(My)-dY; xSearch(Mx)-dX; zSearch(Mz)-dZ] - [lookupTable.yPre+1; lookupTable.xPre+1; lookupTable.zPre{DMD_ix}+1];
    
        % motion = shiftsCenter' + [shifts(rr)-dR shifts(cc)-dC];
    % else %the optimum is at an edge of search range; no superresolution
    %     motionDS(DSframeIx,:) = [ySearch(My); xSearch(Mx); zSearch(Mz)] - [lookupTable.yPre+1; lookupTable.xPre+1; lookupTable.zPre{DMD_ix}+1];
    % end

    brightnessDS(DSframeIx) = scalingFactorTable(My, Mx, Mz);
    % dataMatrix(:,DSframeIx) = data;
    % expectedMatrix(:,DSframeIx) = brightnessDS(DSframeIx) .* lookupTable.likelihood_means{DMD_ix}(ySearch(My), xSearch(Mx), zSearch(Mz),:);
    
    if params.saveTiffs
        for zIdx = 1:numFastZs
            A1 = nan(dmdPixelsPerColumn,dmdPixelsPerRow);
            for cIdx = unique(spCols)'
                spIdxs = find(spCols == cIdx & spPlanes == zIdx);
                % queriedRows = spRows(spIdxs)'+motionDS(DSframeIx,1);
                % rowSpacings = diff(queriedRows);
                
                if numel(spIdxs) > 1
                    A1(:,cIdx+round(motionDS(DSframeIx,2))) = interp1(spRows(spIdxs)'+round(motionDS(DSframeIx,1)),data(spIdxs,1)./spCt(spIdxs).*100,1:dmdPixelsPerColumn);
                else
                    continue;
                end
            end
            fTIF.WriteIMG(uint16(A1));
            if numChannels==2
                A2 = nan(dmdPixelsPerColumn,dmdPixelsPerRow);
                for cIdx = unique(spCols)'
                    spIdxs = find(spCols == cIdx & spPlanes == zIdx);
                    % queriedRows = spRows(spIdxs)'+motionDS(DSframeIx,1);
                    % rowSpacings = diff(queriedRows);
                    if numel(spIdxs) > 1
                        A2(:,cIdx+round(motionDS(DSframeIx,2))) = interp1(spRows(spIdxs)'+round(motionDS(DSframeIx,1)),data(spIdxs,2)./spCt(spIdxs).*100,1:dmdPixelsPerColumn);
                    else
                        continue;
                    end
                end
                fTIF.WriteIMG(uint16(A2));
            end
        end
    end
    
    ySearch = max(1,round(motionDS(DSframeIx,1)+lookupTable.yPre+1) - searchRadius):min(xMotRange,round(motionDS(DSframeIx,1)+lookupTable.yPre+1) + searchRadius);
    xSearch = max(1,round(motionDS(DSframeIx,2)+lookupTable.xPre+1) - searchRadius):min(yMotRange,round(motionDS(DSframeIx,2)+lookupTable.xPre+1) + searchRadius);
    zSearch = max(1,round(motionDS(DSframeIx,3)+lookupTable.zPre{DMD_ix}+1) - searchRadiusZ):min(zMotRange,round(motionDS(DSframeIx,3)+lookupTable.zPre{DMD_ix}+1) + searchRadiusZ);
end
catch ME
    disp(ME);
    registrationFailed = true;
end

if params.saveTiffs
    fTIF.close;
end


if registrationFailed
    disp(['REGISTRATION ERROR OCCURRED FOR FILE: ' fn newline 'YOU MAY NEED TO QC THIS FILE!' newline 'CONTINUING...'])
    return
end

disp(['aligning done - took ' num2str(toc) ' sec'])

%% Upsample if downsampled
% motion = interp1(dsFrames,motionDS,firstLine:lastLine,'linear','extrap');
% 
% brightness = interp1(dsFrames,brightnessDS,frames,'linear','extrap');

%% Save out data

aData.numChannels = numChannels;
aData.frametime = 1/aData.alignHz;
aData.DSframes = DSframes;
aData.motionDSc = motionDS(:,2);
aData.motionDSr = motionDS(:,1);
aData.motionDSz = motionDS(:,3);
aData.brightnessDS = brightnessDS;
aData.logLikelihoodDS = loglikelihoodDS;

disp('Getting online motion correction offsets')
[aData.onlineXshift, aData.onlineYshift, aData.onlineZshift] = getOnlineMotion(hLowLevelDataFile, DSframes);

save(fnAdata, 'aData');

if params.saveTiffs && params.efficientTiffSave
    copyfile(fnwriteTmp,fnwrite);
    delete(fnwriteTmp);
end

end

function meta = loadMetadata(datFilename)
    ix = strfind(datFilename, 'DMD'+digitsPattern(1));
    metaFilename = [datFilename(1:ix+3) '.meta'];
    meta = load(metaFilename, '-mat');
end