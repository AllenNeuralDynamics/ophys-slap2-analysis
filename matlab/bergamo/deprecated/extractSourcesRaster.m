function extractSourcesRaster
dr  ='Z:\scratch\ophys\Adrian\automatic segmentation\708303\11-29-23\scans';

[fn dr] = uigetfile([dr filesep '*REGISTERED*.tif']);

%load high-res image
A = ScanImageTiffReader([dr filesep fn]);
IM = A.data;

%load alignment data
fnStemEnd = strfind(fn, '_REGISTERED') -1;
load([dr filesep fn(1:fnStemEnd) '_ALIGNMENTDATA.mat'], 'aData');

nCh = aData.numChannels;
%reshape channels
IM = reshape(IM, size(IM,1), size(IM,2), nCh, size(IM,3)./nCh);

IMf = squeeze(IM(:,:,1,:));

%discard motion frames
% tmp = aData.aError-smoothdata(aData.aError,2, 'movmedian', ceil(2/aData.frametime));
% discardFrames = tmp>(3*std(tmp));
%IMf(:,:,discardFrames) = nan;

[IMc, peaks] = localizeFlashes(IMf);

%generate an initialization for NMF!
upsample = 8;
rGrid = linspace(1,size(IMc,1), upsample*size(IMc,1)+1);
cGrid = linspace(1,size(IMc,2), upsample*size(IMc,2)+1);

[RR,CC] = ndgrid(rGrid,cGrid);

end