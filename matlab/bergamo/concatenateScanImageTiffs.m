function concatenateScanImageTiffs

[selectFiles, dr] = uigetfile("*.tif","MultiSelect","on");
%%
frameStream = {};
frameDesc = {};

disp([num2str(length(selectFiles)) ' files to concatenate']);

for i = 1:length(selectFiles)
    filename = fullfile(dr,sprintf([selectFiles{1}(1:end-9) '%05d.tif'],i));

    A = ScanImageTiffReader(filename);
    data = A.data;
    for idx = 1:size(data,3)
        frameStream{end+1} = data(:,:,idx);
    end

    frameDesc = [frameDesc; A.descriptions()];
    
    if i == 1
        meta = A.metadata;
        [startIdx, ~] = regexp(meta, 'SI.hScan2D.logFramesPerFile');
        b4 = meta(1:startIdx-1);
        endIdx = regexp(meta, 'SI.hScan2D.logFramesPerFileLock');
        after = meta(endIdx:end);
        
        newMeta = [b4 sprintf('SI.hScan2D.logFramesPerFile = Inf\n') after];
    end
end

hdrSplit = strsplit(newMeta, '\n\n');
    
hdrBuf = [uint8(hdrSplit{1}) 0];
hdrBufLen = length(hdrBuf);

hdrBuf = [hdrBuf uint8([sprintf('\n') hdrSplit{2}]) 0];

pfix = [1 3 3 7 typecast(uint32(4),'uint8') typecast(uint32(hdrBufLen),'uint8') typecast(uint32(length(hdrBuf)-hdrBufLen),'uint8')];

tifHeaderData = [pfix hdrBuf]';
tifHeaderStringOffset = length(pfix);
tifRoiDataStringOffset = length(pfix) + hdrBufLen;
%%
[header,~,~,~] = scanimage.util.opentif(fullfile(dr,sprintf([selectFiles{1}(1:end-9) '%05d.tif'],1)));
[hMroiRoiGroup,~,~] = scanimage.util.readTiffRoiData(fullfile(dr,sprintf([selectFiles{1}(1:end-9) '%05d.tif'],1)));

dataSigned    = true;
bitsPerSample = 16;

numChannelSave = numel(header.SI.hChannels.channelSave);
pixelsPerLine = header.SI.hRoiManager.pixelsPerLine;
linesPerFrame = header.SI.hRoiManager.linesPerFrame;

blankFrameDescription = repmat(' ',1,2000);

imageSize = pixelsPerLine * linesPerFrame * (bitsPerSample/8);

flybackPeriods = ceil(header.SI.hScan2D.flybackTimePerFrame * header.SI.hScan2D.scannerFrequency);
flybackLinesPerFrame = round((flybackPeriods * 2^header.SI.hScan2D.bidirectional)/2)*2;

rois = hMroiRoiGroup.rois;

zs = header.SI.hStackManager.zs;

scanFields = arrayfun(@(z)hMroiRoiGroup.scanFieldsAtZ(z),...
                zs,'UniformOutput',false);
            
cumPixelResolutionAtZ = zeros(0,2);
mRoiLogging = false;
for zidx = 1:length(scanFields)
    sfs = scanFields{zidx};
    pxRes = zeros(0,2);
    for sfidx = 1:length(sfs)
        sf = sfs{sfidx};
        pxRes(end+1,:) = sf.pixelResolution(:)';
    end
    mRoiLogging = mRoiLogging || size(pxRes,1) > 1;
    cumPixelResolutionAtZ(end+1,:) = [max(pxRes(:,1)),sum(pxRes(:,2))+((size(pxRes,1)-1)*flybackLinesPerFrame)];
end

mRoiLogging = mRoiLogging || any(cumPixelResolutionAtZ(1,1) ~= cumPixelResolutionAtZ(:,1));
mRoiLogging = mRoiLogging || any(cumPixelResolutionAtZ(1,2) ~= cumPixelResolutionAtZ(:,2));
linesPerFrame = max(cumPixelResolutionAtZ(:,2));
pixelsPerLine = max(cumPixelResolutionAtZ(:,1));

sf = scanFields{1}{1};
resDenoms = 2^30 ./ (1e4 * sf.pixelResolutionXY ./ (sf.sizeXY * header.SI.objectiveResolution));

xResolutionNumerator = 2^30;
xResolutionDenominator = resDenoms(1);
yResolutionNumerator = 2^30;
yResolutionDenominator = resDenoms(2);
            
%%
hNewTif = scanimage.components.scan2d.TiffStream;

assert(hNewTif.open(fullfile(dr,sprintf([selectFiles{1}(1:end-10) '.tif'])),tifHeaderData,tifHeaderStringOffset,tifRoiDataStringOffset), 'Failed to create log file.');

hNewTif.configureImage(pixelsPerLine, linesPerFrame, (bitsPerSample/8), numChannelSave, dataSigned, blankFrameDescription,...
    xResolutionNumerator, xResolutionDenominator, yResolutionNumerator, yResolutionDenominator);

bitsPerSample = 16;
imageSize = size(data,1) * size(data,2) * (bitsPerSample/8);

for idx = 1:numel(frameStream)
    hNewTif.replaceImageDescription(frameDesc{idx});
    hNewTif.appendFrame(int16(frameStream{idx}), imageSize);
end

hNewTif.close();
hNewTif.cleanUp();  