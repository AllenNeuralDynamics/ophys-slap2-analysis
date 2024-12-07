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
fnRegDS = cell(nDMDs,length(trialTable.trueTrialIx));
fnAdata = cell(nDMDs,length(trialTable.trueTrialIx));
parfor p_ix = 1:numel(fixs)
    f_ix = fixs(p_ix); DMD_ix = dixs(p_ix);
    [fnRegDS{p_ix}, fnAdata{p_ix}]= alignIntegrationAsync(dr, trialTable, lookupTable, params, f_ix, DMD_ix);
end
trialTable.fnRegDSInt = fnRegDS;
trialTable.fnAdataInt = fnAdata;
trialTable.alignParamsInt = params;
save([dr filesep fn], "trialTable")

disp('done integrationRegistration.')
end

function lookupTable = makeRefLookupTable(dr, trialTable, params)
nDMDs = size(trialTable.filename,1);
for DMDix = nDMDs:-1:1
    [~,n] = fileparts(trialTable.filename{DMDix,1});
    n_base = regexprep(n,'-TRIAL[0-9]+$','','ignorecase');
    metaDataFileName = fullfile(dr, [n_base '.meta']);
    mustBeFile(metaDataFileName);
    metaData = load(metaDataFileName, '-mat');

    list = dir([dr filesep '**' filesep '*DMD' int2str(DMDix) '_CONFIG2-REFERENCE*']);
    ReferenceStack_ = slap2.gui.refstack.ReferenceStack.loadTif(fullfile(list(1).folder, list(1).name));
    numChannelsRefStack = numel(ReferenceStack_.data);
    for chIx = numChannelsRefStack:-1:1
        refStack(:,:,:,chIx) = permute(ReferenceStack_.data{chIx},[2 1 3]);
    end

    numLinesPerCycle = length(metaData.AcquisitionContainer.ParsePlan.acqParsePlan);

    dmdPixelsPerColumn = metaData.dmdPixelsPerColumn;
    dmdPixelsPerRow = metaData.dmdPixelsPerRow;
    numFastZs = length(metaData.AcquisitionContainer.ParsePlan.zs);

    % get list of superpixels and extract data

    allSuperPixelIDs{DMDix} = [];

    fastZ2RefZ = zeros(numFastZs,1);
    for z = 1:numFastZs
        [~, ind] = min(abs(ReferenceStack_.zs - metaData.AcquisitionContainer.ParsePlan.zs(z)));
        fastZ2RefZ(z) = ind;
    end

    for lineSweepIdx = 1:numLinesPerCycle
        superPixIdxs = metaData.AcquisitionContainer.ParsePlan.acqParsePlan(lineSweepIdx).superPixelID;

        if size(superPixIdxs,1) == 0; continue; end

        zIdx = metaData.AcquisitionContainer.ParsePlan.acqParsePlan(lineSweepIdx).sliceIdx(1) + 1;

        % only take integration mode pixels
        superPixIdxs(superPixIdxs <= dmdPixelsPerColumn*dmdPixelsPerRow) = [];

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
    allPixelReplacementMaps = metaData.AcquisitionContainer.ParsePlan.pixelReplacementMaps;

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
        floor((max(size(bl,3),numFastZs)-size(bl,3))/2)],...
        mean(bl(:)),...
        'both');

    if size(bl,1) < dmdPixelsPerColumn
        bl = padarray(bl,[1,0,0],mean(bl(:)),'pre');
    end

    if size(bl,2) < dmdPixelsPerRow
        bl = padarray(bl,[0,1,0],mean(bl(:)),'pre');
    end

    % Calculate lookup table

    xPre = params.maxshiftXY; xPost = params.maxshiftXY;
    yPre = params.maxshiftXY; yPost = params.maxshiftXY;
    zPre = params.maxshiftZ; zPost = params.maxshiftZ;

    likelihood_means{DMDix} = makeLookupTable(bl, sparseMaskInds{DMDix}, numFastZs, fastZ2RefZ,-yPre:yPost,-xPre:xPost,-zPre:zPost,ReferenceStack_.channels);
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
hLowLevelDataFile = hSlap2DataFile.hDataFile;
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
    disp(['Requested DS freq too fast, adjusting to ' aData.alignHz ' Hz']);
end

fnwrite = [dr filesep fnW '_REGISTERED_DOWNSAMPLED-' int2str(aData.alignHz) 'Hz.tif'];
if params.efficientTiffSave
    fnwriteTmp = [params.tempFileDir filesep int2str(round(rand(1)*10000)) '.tif'];
else
    fnwriteTmp = fnwrite;
end

fnAdata = [dr filesep fnW '_ALIGNMENTDATA.mat'];

if ~params.overwriteExisting && exist(fnAdata, 'file') && exist(fnwrite, 'file')
    disp([fn ' is already aligned; skipping' newline 'To force realign, pass TRUE as second argument']);
    return
end

DSframes = ceil(firstLine:dt:lastLine);
nDSframes= length(DSframes); %number of downsampled frames

%% load metadata
dmdPixelsPerColumn = hLowLevelDataFile.metaData.dmdPixelsPerColumn;
dmdPixelsPerRow = hLowLevelDataFile.metaData.dmdPixelsPerRow;
numFastZs = length(hLowLevelDataFile.fastZs);

%%  Calculate log likelihoods and infer motion from MAP
log_means = log(lookupTable.likelihood_means{DMD_ix} + 1e-8);

% how much change is allowed in each dimension at each step
searchRadius = aData.clipShift;

motionDS = nan(nDSframes,3);
brightnessDS = nan(nDSframes,1);
loglikelihoodDS = nan(nDSframes,1);

% dataMatrix = zeros(length(lookupTable.allSuperPixelIDs{DMD_ix}),nDSframes);
% expectedMatrix = zeros(length(lookupTable.allSuperPixelIDs{DMD_ix}),nDSframes);

pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM; 250nm
fTIF = Fast_BigTiff_Write(fnwriteTmp,pixelscale,0);
    
xMotRange = lookupTable.xPre + lookupTable.xPost + 1;
yMotRange = lookupTable.yPre + lookupTable.yPost + 1;
zMotRange = lookupTable.zPre + lookupTable.zPost + 1;

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
[spRows, spCols] = ind2sub([dmdPixelsPerColumn, dmdPixelsPerRow], medianIndices);

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
    allLineData = hLowLevelDataFile.getLineData(lineIndices, cycleIndices);
    for t = 1:length(allLineData)
        superPixIdxs = hLowLevelDataFile.lineSuperPixelIDs{lineIndices(t)};

        if size(superPixIdxs,1) == 0; continue; end

        lineData = allLineData{t};
        zIdx = hLowLevelDataFile.lineFastZIdxs(lineIndices(t));

        spID = superPixIdxs*100 + uint32(zIdx); % make superpixel index with Z plane
        [~,spIdx] = ismember(spID,lookupTable.allSuperPixelIDs{DMD_ix});

        data(spIdx(spIdx>0),:) = data(spIdx(spIdx>0),:) + single(lineData(spIdx>0,:));
        spCt(spIdx(spIdx>0)) = spCt(spIdx(spIdx>0)) + 1;
    end
    data = data ./ 100; %./ spSampleCt;
    data(spCt == 0,:) = nan;

    if mean(spCt == 0) > 0.5; continue; end

    % calculate log likelihoods at all shifts
    [logLikelihoodTable, scalingFactorTable] = poissonLogLikelihoodTable(data, lookupTable.likelihood_means{DMD_ix} ...
                                                                                .* repmat(reshape(spCt,[1 1 1 1 length(lookupTable.allSuperPixelIDs{DMD_ix})]),[size(lookupTable.likelihood_means{DMD_ix},1:4) 1]), ...
                                                                        log_means,ySearch,xSearch,zSearch,channels,params.robust);

    [loglikelihoodDS(DSframeIx), I] = max(logLikelihoodTable(:));

    [My, Mx, Mz] = ind2sub(size(logLikelihoodTable),I);

    if My>1 && My<size(logLikelihoodTable,1) && Mx>1 && Mx<size(logLikelihoodTable,2) && Mz>1 && Mz<size(logLikelihoodTable,3)
        %perform superresolution upsampling, assuming quadratic loglikelihood
        ratioY = min(1e6,(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My-1,Mx,Mz))/(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My+1,Mx,Mz)));
        dY = (1-ratioY)/(1+ratioY)/2;
        
        ratioX =min(1e6, (logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx-1,Mz))/(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx+1,Mz)));
        dX = (1-ratioX)/(1+ratioX)/2;

        ratioZ =min(1e6, (logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx,Mz-1))/(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx,Mz+1)));
        dZ = (1-ratioZ)/(1+ratioZ)/2;

        motionDS(DSframeIx,:) = [ySearch(My)-dY; xSearch(Mx)-dX; zSearch(Mz)-dZ] - [lookupTable.yPre+1; lookupTable.xPre+1; lookupTable.zPre+1];
    
        % motion = shiftsCenter' + [shifts(rr)-dR shifts(cc)-dC];
    else %the optimum is at an edge of search range; no superresolution
        motionDS(DSframeIx,:) = [ySearch(My); xSearch(Mx); zSearch(Mz)] - [lookupTable.yPre+1; lookupTable.xPre+1; lookupTable.zPre+1];
    end

    brightnessDS(DSframeIx) = scalingFactorTable(My, Mx, Mz);
    % dataMatrix(:,DSframeIx) = data;
    % expectedMatrix(:,DSframeIx) = brightnessDS(DSframeIx) .* lookupTable.likelihood_means{DMD_ix}(ySearch(My), xSearch(Mx), zSearch(Mz),:);
    
    A1 = nan(dmdPixelsPerColumn,dmdPixelsPerRow);
    for cIdx = unique(spCols)'
        spIdxs = find(spCols == cIdx);
        % queriedRows = spRows(spIdxs)'+motionDS(DSframeIx,1);
        % rowSpacings = diff(queriedRows);
        
        if numel(spIdxs) > 1
            A1(:,cIdx+round(motionDS(DSframeIx,2))) = interp1(spRows(spIdxs)'+round(motionDS(DSframeIx,1)),data(spIdxs,1)./spCt(spIdxs),1:dmdPixelsPerColumn);
        else
            continue;
        end
    end
    fTIF.WriteIMG(single(A1));
    if numChannels==2
        A2 = nan(dmdPixelsPerColumn,dmdPixelsPerRow);
        for cIdx = unique(spCols)'
            spIdxs = find(spCols == cIdx);
            % queriedRows = spRows(spIdxs)'+motionDS(DSframeIx,1);
            % rowSpacings = diff(queriedRows);
            if numel(spIdxs) > 1
                A2(:,cIdx+round(motionDS(DSframeIx,2))) = interp1(spRows(spIdxs)'+round(motionDS(DSframeIx,1)),data(spIdxs,2)./spCt(spIdxs),1:dmdPixelsPerColumn);
            else
                continue;
            end
        end
        fTIF.WriteIMG(single(A2));
    end
    
    ySearch = max(1,round(motionDS(DSframeIx,1)+lookupTable.yPre-1) - searchRadius):min(xMotRange,round(motionDS(DSframeIx,1)+lookupTable.yPre-1) + searchRadius);
    xSearch = max(1,round(motionDS(DSframeIx,2)+lookupTable.xPre-1) - searchRadius):min(yMotRange,round(motionDS(DSframeIx,2)+lookupTable.xPre-1) + searchRadius);
    zSearch = max(1,round(motionDS(DSframeIx,3)+lookupTable.zPre-1) - searchRadius):min(zMotRange,round(motionDS(DSframeIx,3)+lookupTable.zPre-1) + searchRadius);
end
catch ME
    disp(ME);
    registrationFailed = true;
end

fTIF.close;


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

if params.efficientTiffSave
    copyfile(fnwriteTmp,fnwrite);
    delete(fnwriteTmp);
end

end

function meta = loadMetadata(datFilename)
    ix = strfind(datFilename, 'DMD'+digitsPattern(1));
    metaFilename = [datFilename(1:ix+3) '.meta'];
    meta = load(metaFilename, '-mat');
end