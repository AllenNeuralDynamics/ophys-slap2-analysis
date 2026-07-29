function alignSingleFileScanImageStack(dr, fn, options)

arguments
    dr = []
    fn = []
    options.pixelSizeUm (1,1) double = 0
    options.sqrtTform (1,1) logical = false
end

pixelSizeUm = options.pixelSizeUm;
sqrtTform = options.sqrtTform;

if isempty(dr) || isempty(fn)
    [fn, dr] = uigetfile('*.tif*', ...
        'Select the individual full stack file', ...
        'multiselect', 'off');
end

tiffFn = fullfile(dr,fn);
[IMs,meta] = readTiff(tiffFn);

metaLines = strsplit(meta, '\n');
for lineIx = 1:length(metaLines)
    try
        eval([metaLines{lineIx} ';']);
    catch ME
        continue
    end
end

numChannels = length(SI.hChannels.channelSave);
numFramesPerSlice = SI.hStackManager.framesPerSlice;
zStep = SI.hStackManager.actualStackZStepSize;
numZs = SI.hStackManager.actualNumSlices;

if pixelSizeUm <= 0
    try
        pixelSizeUm = abs(diff(SI.hRoiManager.imagingFovUm(1:2,1)))./SI.hRoiManager.pixelsPerLine;
    catch
        disp('Error calculating pixel size, using 0.5um');
        pixelSizeUm = 0.5;
    end
end

[b,a] = butter(4, 0.05, 'high'); %butterworth filter for computing correlations
IM = []; IMc = []; IMsk = [];
for Zix = numZs:-1:1
    disp(['Aligning plane: ' int2str(Zix) ' of ' int2str(numZs)])
    imIdxsThisPlane = numChannels*numFramesPerSlice*(Zix-1)+ (1:(numChannels*numFramesPerSlice));
    
    % data order is XYCT
    IMtmp = permute(single(IMs(:,:,imIdxsThisPlane)), [2 1 3]); %transpose image; needed for bidi artifact detection by normcorr
    IMtmp = reshape(IMtmp,size(IMtmp,1),size(IMtmp,2),numChannels,numFramesPerSlice);

    % align
    [IM(:,:,:,Zix), IMc(:,:,:,Zix), IMsk(:,:,:,Zix)] = alignMultiChannel(IMtmp, b, a);
end

%straighten out the reference stack
IM = straightenStack(IM, IMc, IMsk); 

% %bin and square root
% IM = IM(1:2:end,1:2:end) + IM(2:2:end,1:2:end) +IM(1:2:end,2:2:end) + IM(2:2:end,2:2:end);
% downsampledStack = arrayfun(@(x) IM(1:2:end, 1:2:end, x) + IM(2:2:end, 1:2:end, x) + IM(1:2:end, 2:2:end, x) + IM(2:2:end, 2:2:end, x), 1:size(IM, 3), 'UniformOutput', false);
% downsampledStack = cat(3, downsampledStack{:});

if sqrtTform
    IM = max(0, IM- prctile(IM(:), 1));
    IM = sqrt(max(0,IM));
end

%remap refIM to full dynamic range; off by default
normalizeIntensity = false;
if normalizeIntensity
    IM = IM./max(IM, [], 'all'); %./max(IMout, 'all')
    IM = uint16(IM.*65000); 
else
    IM = uint16(IM);
end

%save reference image, one stack for each channel
disp('Saving...')

outputPathCh1 = [tiffFn(1:end-4) '-REF_Ch1.ome.tif'];
%outputPathCh2 = [tiffFn(1:end-4) '-REF_Ch2.ome.tif'];
%outputPathCh1corr = [tiffFn(1:end-4) '-CORR_Ch1.ome.tif'];
%outputPathCh2corr = [tiffFn(1:end-4) '-CORR_Ch2.ome.tif'];
%outputPathCh1sk = [tiffFn(1:end-4) '-SK_Ch1.ome.tif'];
%outputPathCh2sk = [tiffFn(1:end-4) '-SK_Ch2.ome.tif'];

bfCheckJavaPath;
metadata1 = createMinimalOMEXMLMetadata(squeeze(IM(:,:,1,:)));
%metadata2 = createMinimalOMEXMLMetadata(squeeze(IM(:,:,2,:)));
pixelSizeObj = ome.units.quantity.Length(java.lang.Double(pixelSizeUm), ome.units.UNITS.MICROMETER);
metadata1.setPixelsPhysicalSizeX(pixelSizeObj,0); %metadata2.setPixelsPhysicalSizeX(pixelSizeObj, 0);
metadata1.setPixelsPhysicalSizeY(pixelSizeObj, 0); %metadata2.setPixelsPhysicalSizeY(pixelSizeObj, 0);
pixelSizeZObj = ome.units.quantity.Length(java.lang.Double(abs(zStep)), ome.units.UNITS.MICROMETER);
metadata1.setPixelsPhysicalSizeZ(pixelSizeZObj, 0); %metadata2.setPixelsPhysicalSizeZ(pixelSizeZObj, 0);
bfsave(squeeze(IM(:,:,1,:)), outputPathCh1, 'BigTiff', true, 'metadata', metadata1);
%bfsave(squeeze(IM(:,:,2,:)), outputPathCh2, 'BigTiff', true, 'metadata', metadata2);
%bfsave(squeeze(IMc(:,:,1,:)), outputPathCh1corr, 'BigTiff', true);
%bfsave(squeeze(IMc(:,:,2,:)), outputPathCh2corr, 'BigTiff', true);
%bfsave(squeeze(IMsk(:,:,1,:)), outputPathCh1sk, 'BigTiff', true);
%bfsave(squeeze(IMsk(:,:,2,:)), outputPathCh2sk, 'BigTiff', true);



end



