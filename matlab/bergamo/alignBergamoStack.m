%alignBergamoStack
%Aligns a scanimage slow stack with repeats acquired on the Bergamo, generating a single
%high-quality image
%Stacks should be 1024 x 1024
%!!! The stack should be split into 1 file per Z-position. This reduces memory load


[fn, dr] = uigetfile('*.tif*', 'Select ONE FILE PER PLANE', 'multiselect', 'on');
if isnumeric(fn)
    return %
    % 
    % user abort
end
if ~iscell(fn)
    fn = {fn};
end
tiffFn = fullfile(dr,fn);
tiffFn = sort_nat(tiffFn);

[IMs,meta] = readTiff(tiffFn{1});
metaLines = strsplit(meta, '\n');
for lineIx = 1:length(metaLines)
    try
        eval([metaLines{lineIx} ';']);
    catch ME
        continue
    end
end

% try
%     eval(meta);
% catch ME
%     disp(ME.message);
% end


numChannels = length(SI.hChannels.channelSave);
numFramesPerSlice = SI.hStackManager.framesPerSlice;
zStep = SI.hStackManager.actualStackZStepSize;
numZs = length(fn);
pixelSizeUm = abs(diff(SI.hRoiManager.imagingFovUm(1:2,1)))./SI.hRoiManager.pixelsPerLine;

[b,a] = butter(4, 0.05, 'high'); %butterworth filter for computing correlations

for Zix = numZs:-1:1
    disp(['Aligning plane: ' int2str(Zix) ' of ' int2str(numZs)])
    [IMs,~] = readTiff(tiffFn{Zix});

    % data order is XYCT
    IMs = permute(single(IMs), [2 1 3]); %transpose image; needed for bidi artifact detection by normcorr
    IMs = reshape(IMs,size(IMs,1),size(IMs,2),numChannels,numFramesPerSlice);

    % align
    [refIM(:,:,:,Zix), IMc(:,:,:,Zix)] = alignMultiChannel(IMs, b, a);
end

%straighten out the reference stack
[IM, IMc] = straightenStack(refIM, IMc);

%convert to 16 bit
IM = max(0, IM- prctile(IM(:), 1));
IMc = max(0, IMc- prctile(IMc(:), 1));
IMc = IMc./max(IMc, [], 'all');
IMc = uint16(IMc.*65000);

%remap to full dynamic range; optional
normalizeIntensity = false;
if normalizeIntensity
    IM = IM./max(IM, [], 'all'); %./max(IMout, 'all')
    IM = uint16(IM.*65000); 
else
    IM = uint16(IM);
end

%save reference image, one stack for each channel
disp('Saving...')

outputPathCh1 = [tiffFn{1}(1:end-4) '-REF_Ch1.ome.tif'];
outputPathCh2 = [tiffFn{1}(1:end-4) '-REF_Ch2.ome.tif'];
outputPathCh1corr = [tiffFn{1}(1:end-4) '-CORR_Ch1.ome.tif'];
outputPathCh2corr = [tiffFn{1}(1:end-4) '-CORR_Ch2.ome.tif'];

bfCheckJavaPath;
metadata1 = createMinimalOMEXMLMetadata(squeeze(IM(:,:,1,:)));
pixelSizeObj = ome.units.quantity.Length(java.lang.Double(pixelSizeUm), ome.units.UNITS.MICROMETER);
metadata1.setPixelsPhysicalSizeX(pixelSizeObj,0);
metadata1.setPixelsPhysicalSizeY(pixelSizeObj, 0);
pixelSizeZObj = ome.units.quantity.Length(java.lang.Double(abs(zStep)), ome.units.UNITS.MICROMETER);
metadata1.setPixelsPhysicalSizeZ(pixelSizeZObj, 0); 
bfsave(squeeze(IM(:,:,1,:)), outputPathCh1, 'BigTiff', true, 'metadata', metadata1);
bfsave(squeeze(IMc(:,:,1,:)), outputPathCh1corr, 'BigTiff', true);

if size(IM,3) == 2
    metadata2 = createMinimalOMEXMLMetadata(squeeze(IM(:,:,2,:)));
    metadata2.setPixelsPhysicalSizeX(pixelSizeObj, 0);
    metadata2.setPixelsPhysicalSizeY(pixelSizeObj, 0);
    metadata2.setPixelsPhysicalSizeZ(pixelSizeZObj, 0);
    bfsave(squeeze(IM(:,:,2,:)), outputPathCh2, 'BigTiff', true, 'metadata', metadata2);
    bfsave(squeeze(IMc(:,:,2,:)), outputPathCh2corr, 'BigTiff', true);
end

function [refPlane, corrPlane] = alignMultiChannel(IMs, b,a)
Y = squeeze(makeHighPass(sum(IMs,3))); % data order is XYCTZV
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',100,'us_fac',4,'init_batch',40, 'correct_bidir', true);
refPlane = [];
[~,shifts,~,options, col_shift] = normcorre(Y,options_rigid);
for chIx= size(IMs,3):-1:1
    tmp= apply_shifts(squeeze(IMs(:,:,chIx,:)),shifts,options,0,0,0,col_shift);
   
    refPlane(:,:,chIx) = mean(tmp,3);

    %compute correlation image
    %highpass in time
    HP = permute(filtfilt(b,a,permute(double(tmp), [3 1 2])), [2 3 1]);
    ss = sum(HP.^2,3, 'omitnan');
    vertC = sum(HP .* circshift(HP, [1 0 0 0]),3, 'omitnan')./sqrt(ss.*circshift(ss, [1 0 0 ]));
    horzC = sum(HP .* circshift(HP, [0 1 0 0]),3, 'omitnan')./sqrt(ss.*circshift(ss, [0 1 0 ]));
    C = mean(cat(3, horzC, circshift(horzC,1,2), vertC, circshift(vertC, 1,1)),3, 'omitnan');
% 
%     %compute skewness image
%     HP = imgaussfilt(HP, [0.6 0.6]);
%     sk = skewness(HP,0,3).*sqrt(max(0, refPlane(:,:,chIx));
%     obj.IMsk = (sk-median(sk(:), 'omitnan'))./std(sk, 0,'all','omitnan');
%     obj.IMact = obj.IMc + obj.IMsk;

    
    corrPlane(:,:,chIx) = C;
end
end



function IMHP = makeHighPass(IM1)
IMHP = IM1; %sum the channels to speed alignment
IMHP(IMHP<0) = 0;
IMHP(isnan(IMHP)) = 0;
IMHP = imgaussfilt(IMHP, [0.5 0.5])-imgaussfilt(IMHP, [4 4]);
end

function [IMs,meta, ImDescriptions] = readTiff(tiffFn)
hTiff= ScanImageTiffReader(tiffFn);
try
    ImDescriptions = hTiff.descriptions;
    meta = hTiff.metadata();
    IMs = hTiff.data;
catch ME
    most.idioms.safeDeleteObj(hTiff);
    ME.rethrow();
end
most.idioms.safeDeleteObj(hTiff);
end
