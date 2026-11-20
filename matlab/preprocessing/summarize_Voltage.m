function summarize_Voltage(dr_or_pathToTrialTable, paramsIn)

import ScanImageTiffReader.ScanImageTiffReader
if nargin>1
    if ischar(paramsIn)  % Parse JSON String to Structure
        paramsIn = jsondecode(paramsIn);
    end
    params = setParams('summarize_Voltage', paramsIn);
else
    params = setParams('summarize_Voltage');
end
if ~nargin
    [trialTablefn, dr] = uigetfile({'*trialTable*.mat;*trial_table*.h5', 'Trial table (*.mat, *.h5)'}, 'Select trial table');
else
    %parse dr
    %_or_pathToTrialTable
    if exist(dr_or_pathToTrialTable, 'dir')
        dr = dr_or_pathToTrialTable;
        if exist(fullfile(dr, 'trial_table.h5'), 'file')
            trialTablefn = 'trial_table.h5'; % GIAnT-MATLAB layout
        else
            trialTablefn = 'trialTable.mat';
        end
    else
        [dr, trialTablefn, ext] = fileparts(dr_or_pathToTrialTable);
        trialTablefn = [trialTablefn ext];
    end
end
copyReadDeleteScanImageTiff([]); %make sure we can use the function in parallel loops

%confirm that all files exist. dirs.dataDr holds raw .dat files and reference
%TIFFs, dirs.mocoDr the alignment data, dirs.saveDr the output root. For a
%trialTable.mat all three equal dr; for a GIAnT trial_table.h5 they come from
%the table's datadr/savedr fields.
[trialTable, keepTrials, dirs] = verifyFiles(trialTablefn, dr, params);
if isfield(trialTable, 'refStack') && iscell(trialTable.refStack)
    for dmdIx = 1:numel(trialTable.refStack)
        if isstruct(trialTable.refStack{dmdIx}) && isfield(trialTable.refStack{dmdIx}, 'IM')
            trialTable.refStack{dmdIx}.IM = []; %this uses a lot of memory and we won't need it
        end
    end
end
nTrials = size(trialTable.filename,2);

disp(['## SUMMARIZING' newline 'Data folder:'])
disp(dirs.dataDr)

savedr = [dirs.saveDr filesep 'ExperimentSummary'];
% on CodeOcean /data is read-only and we save to /results
is_CodeOcean = ~(getenv("CO_CPUS") == "");
if is_CodeOcean
    savedr = strrep(savedr, '/data', '/results');
end
if ~exist(savedr, 'dir')
    mkdir(savedr);
end
fnsave = [savedr filesep 'Summary_Voltage-' datestr(now, 'YYmmDD-HHMMSS') '.mat'];

%trace any manual ROIs
nDMDs =size(trialTable.filename,1);
rTot = 0;
for DMDix = 1:nDMDs
    %load image data
    firstValidTrial = find(~cellfun(@isempty, trialTable.filename(DMDix,:)), 1, 'first');
    hDataFile = slap2.util.DataFile([dirs.dataDr filesep trialTable.filename{DMDix,firstValidTrial}]);

    imagingRois = hDataFile.metaData.AcquisitionContainer.ROIs.rois;
    isIntegration = false(1,length(imagingRois));
    for rix = 1:length(imagingRois)
        isIntegration(rix) = strcmpi(imagingRois{rix}.imagingMode, 'Integrate');
    end
    integrationRois = imagingRois(isIntegration);
    nImagingROIs = numel(integrationRois);

    %GET ROI outlines
    outlines = {};
    maskImage{DMDix} = -1.*ones(hDataFile.metaData.dmdPixelsPerColumn, hDataFile.metaData.dmdPixelsPerRow);
    for rix = 1:nImagingROIs
        shape = integrationRois{rix}.shapeData;
        tmp = false(hDataFile.metaData.dmdPixelsPerColumn, hDataFile.metaData.dmdPixelsPerRow);
        tmp(sub2ind(size(tmp), shape(:,1), shape(:,2))) = true;
        summary.masks{DMDix}(:,:,rix) = tmp;
        if DMDix==1
            maskImage{DMDix}(summary.masks{DMDix}(:,:,rix)) = rix;
        else
            maskImage{DMDix}(summary.masks{DMDix}(:,:,rix)) = rix+size(summary.masks{1},3);
        end
        %make an image where each pixel corresponds to the mask ID
        outlines = cat(1, outlines, bwboundaries(tmp));
    end

    %load reference image
    refFn = dir([dirs.dataDr filesep '**' filesep '*DMD' int2str(DMDix) '*REFERENCE*.tif']);
    assert(length(refFn)==1, 'Found multiple reference images in the folder (possibly in nested folders!)');
    A = ScanImageTiffReader([refFn.folder filesep refFn.name]);
    IDs = A.descriptions;
    for im_ix = length(IDs):-1:1
        js = jsondecode(IDs{im_ix});
        z(im_ix) = js.z;
        ch(im_ix) = js.channel;
    end
    nChan = numel(unique(ch));
    Zs = unique(z);
    stack = A.data();
    stack = reshape(stack, size(stack,1), size(stack,2), nChan, []);

    %get imaging Z plane
    metaZ = hDataFile.metaData.AcquisitionContainer.ParsePlan.zs;
    [~, bestZix] = min(abs(Zs-metaZ)); bestZix = bestZix(1);
    %bestZ = Zs(bestZix);

    %The image plane to trace on
    imPlane = permute(stack(:,:,:,bestZix),[2 1 3]);
    summary.refPlane{DMDix} = imPlane;
    clear stack

    if params.manualRois
        hROIs = drawROIs(imPlane, dirs.saveDr, refFn.name, outlines); %manual ROIs are saved next to the other outputs
        waitfor(hROIs.hF);
        nAnalysisROIs(DMDix) = numel(hROIs.roiData);
        summary.masks{DMDix} = false(0);
        maskImage{DMDix}(:,:) = -1;
        for rix = 1:nAnalysisROIs(DMDix)
            rTot = rTot+1;
            summary.masks{DMDix}(:,:,rix) = hROIs.roiData{rix}.mask;
            %make an image where each pixel corresponds to the mask ID
            maskImage{DMDix}(summary.masks{DMDix}(:,:,rix)) = rTot;
        end
    else
        nAnalysisROIs(DMDix) = nImagingROIs;
    end

    figure('name', ['Reference Image for DMD' int2str(DMDix)]),
    imshow(summary.refPlane{DMDix}(:,:,1),[]);
    figure('name', ['Masks for DMD' int2str(DMDix)]),
    imshow(maskImage{DMDix},[]); colormap('jet')
    drawnow
end
summary.userROIs = maskImage;

for DMDix = 1:nDMDs

    %extract signal for each superpixel
    disp('Reading Voltage data')
    % mIM = cell(1, nTrials);
    % fnR = trialTable.fnRegDS(DMDix, :);
    fns = trialTable.fnRaw(DMDix,:);
    fnA = trialTable.fnAdata(DMDix,:);
    fls = trialTable.firstLine(DMDix,:);
    els = trialTable.lastLine(DMDix,:);
    %lookup = makeOrLoadLookupTable(dr, trialTable, DMDix, params);

    E = cell(1, nTrials);
    aData = cell(1, nTrials);
    masks = summary.masks{DMDix};
    summary.footprints{DMDix} = [];
    for trialIx = 1:nTrials
        disp(['Processing Trial ' int2str(trialIx)]);
        if keepTrials(DMDix,trialIx)
            aData{trialIx} = readAlignmentData(fullfile(dirs.mocoDr, fnA{trialIx})); %LOAD ALIGNMENT DATA (.mat or GIAnT .h5)
            if isfield(aData{trialIx}, 'varFacDS')
                aData{trialIx} = rmfield(aData{trialIx}, 'varFacDS');
            end
            [~,fn, ext] = fileparts(fns{trialIx});
            if isempty(summary.footprints{DMDix})
                [E{trialIx}, summary.footprints{DMDix}]= loadAndProcessVoltage(dirs.dataDr, [fn ext], fls(trialIx), els(trialIx), aData{trialIx}, masks, params);
            else
                E{trialIx}= loadAndProcessVoltage(dirs.dataDr, [fn ext], fls(trialIx), els(trialIx), aData{trialIx}, masks, params);
            end
        end
    end
    summary.E(:,DMDix) = E;
    summary.aData{:,DMDix} = aData;
end

%prepare file for saving
summary.params = params;
summary.trialTable = trialTable;
summary.dr = dirs.dataDr;
summary.dirs = dirs;
summary.trialTablePath = fullfile(dr, trialTablefn);

%save
save(fnsave, 'summary', "-v7.3");
disp('Done summarize_Voltage')
end



% function lookupTable = makeOrLoadLookupTable(dr, trialTable, DMDix, params)
% [~,n] = fileparts(trialTable.filename{DMDix,1});
% n_base = regexprep(n,'-TRIAL[0-9]+$','','ignorecase');
% metaDataFileName = fullfile(dr, [n_base '.meta']);
% mustBeFile(metaDataFileName);
% metaData = load(metaDataFileName, '-mat');
% 
% if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
%     numLinesPerCycle = length(metaData.AcquisitionContainer.AcquisitionPlan.superPixelIDs);
% 
%     zs_ix = horzcat(metaData.AcquisitionContainer.AcquisitionPlan.activeZs{:});
%     zs_ix = unique(zs_ix);
% 
%     zPlanes = nan(1,numLinesPerCycle);
%     for ix = 1:numLinesPerCycle
%         if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan.activeZs{ix})
%             zPlanes(ix) = metaData.AcquisitionContainer.AcquisitionPlan.activeZs{ix}(1);
%         end
%     end
%     for ix = 1:length(zs_ix)
%         zPlane_um = mean(metaData.AcquisitionContainer.AcquisitionPlan.zTrajectory(zPlanes == zs_ix(ix)));
%         zs(ix) = zPlane_um;
%     end
% else
%     zs = metaData.AcquisitionContainer.ParsePlan.zs;
%     numLinesPerCycle = length(metaData.AcquisitionContainer.ParsePlan.acqParsePlan);
% end
% 
% dmdPixelsPerColumn = metaData.dmdPixelsPerColumn;
% dmdPixelsPerRow = metaData.dmdPixelsPerRow;
% numFastZs = length(zs);
% 
% % get list of superpixels and extract data
% allSuperPixelIDs = [];
% 
% for lineSweepIdx = 1:numLinesPerCycle
%     if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
%         superPixIdxs = metaData.AcquisitionContainer.AcquisitionPlan.superPixelIDs{lineSweepIdx}';
%     else
%         superPixIdxs = metaData.AcquisitionContainer.ParsePlan.acqParsePlan(lineSweepIdx).superPixelID;
%     end
% 
%     if numel(superPixIdxs) == 0; continue; end
% 
%     if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
%         zIdx = metaData.AcquisitionContainer.AcquisitionPlan.activeZs{lineSweepIdx}(1);
%     else
%         zIdx = metaData.AcquisitionContainer.ParsePlan.acqParsePlan(lineSweepIdx).sliceIdx(1)+1;
%     end
% 
%     % only take integration mode pixels
%     if params.integrationOnly
%         superPixIdxs(superPixIdxs <= dmdPixelsPerColumn*dmdPixelsPerRow) = [];
%     end
% 
%     spIDs = superPixIdxs*100+zIdx; % add z plane to end of superpixel ID
%     allSuperPixelIDs = [allSuperPixelIDs; spIDs]; % make list of unique superpixels across all Zs
% end
% 
% [allSuperPixelIDs, ~, ic] = unique(allSuperPixelIDs);
% spSampleCt = accumarray(ic, 1);
% 
% fprintf("%d superpixels detected\n", length(allSuperPixelIDs));
% 
% % make sparse matrix with each superpixel's corresponding mask (roiMasks)
% fprintf("Calculating ROI Masks... ");
% tic;
% 
% % using sparse matrix
% sparseMaskInds = [];
% if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
%     allPixelReplacementMaps = metaData.AcquisitionContainer.AcquisitionPlan.pixelReplacementMaps;
% else
%     allPixelReplacementMaps = metaData.AcquisitionContainer.ParsePlan.pixelReplacementMaps;
% end
% 
% for i = 1:length(allSuperPixelIDs)
%     %tmpMask = zeros(dmdPixelsPerColumn,dmdPixelsPerRow,numFastZs);
% 
%     sp = allSuperPixelIDs(i);
%     zIdx = rem(sp,100);
%     superPixIdx = (sp - zIdx) / 100;
%     pixelReplacementMap = allPixelReplacementMaps{zIdx};
% 
%     open = uint32(pixelReplacementMap(pixelReplacementMap(:,2) == superPixIdx,1))+1;
% 
%     if isempty(open)
%         open = superPixIdx+1;
%     end
% 
%     openR = idivide(open-1, dmdPixelsPerRow, 'floor')+1;
%     openC = open - (openR-1) * dmdPixelsPerRow;
%     openPixs = uint32(openR + (openC-1) * dmdPixelsPerColumn + double(zIdx-1) * dmdPixelsPerColumn * dmdPixelsPerRow);
% 
%     sparseMaskInds = [sparseMaskInds; openPixs, ones(size(openPixs))*i];
% end
% %clear('tmpMask');
% roiMasks = sparse(sparseMaskInds(:,1),sparseMaskInds(:,2),1,dmdPixelsPerColumn * dmdPixelsPerRow * numFastZs,length(allSuperPixelIDs));
% fprintf('done LUT- took %f sec\n',  toc);
% 
% lookupTable.allSuperPixelIDs = allSuperPixelIDs;
% lookupTable.spSampleCt = spSampleCt;
% %lookupTable.sparseMaskInds = sparseMaskInds;
% lookupTable.roiMasks = roiMasks;
% end

