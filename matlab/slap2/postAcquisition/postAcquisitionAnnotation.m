function postAcquisitionAnnotation (dr)
%Runs the soma annotater on raw acquisition files and saves resulting annotations;
% these have to be aligned to the later concensus in summarize_Loco to be used
if ~nargin
    dr = pwd;
end
import ScanImageTiffReader.ScanImageTiffReader;

nDMDs = 2;
fn{1} = uigetfile('*DMD1*.tif', 'Select a DMD1 acquisition file for ANNOTATION');
fn{2} = uigetfile('*DMD2*.tif', 'Select a DMD2 acquisition file for ANNOTATION');

fnAnn = [dr filesep 'postAcqANNOTATIONS.mat'];

if exist(fnAnn, 'file')
    abort = ~strcmpi(questdlg('Post-acquisition annotations already exist. Continue?'), 'Yes');
    if abort
        return
    end
end

for DMDix = 1:nDMDs
    %load image data
    A = ScanImageTiffReader(fullfile(dr,fn{DMDix}));
    IM  = A.data;
    desc = A.descriptions();
    meta = jsondecode(desc{1});
    nchan = meta.numChannels;
    nframes = size(IM,3)./nchan;

    nans = repmat(all((IM==0),3), 1,1,nchan); %a hack to quickly guess which pixels are unimaged
    IM = double(reshape(IM, size(IM,1), size(IM,2), nchan, nframes));
    
    IMt = squeeze(trimmean(IM, 40, 4));
    IMt(nans) = nan;
    IMt = permute(IMt, [2 1 3]);

    hROIs(DMDix) = drawROIs(sqrt(max(0,IMt(:,:,end:-1:1))), dr, fn{DMDix}); %3rd dim inverted because channels are GR and colors are RGB

    ROIs(DMDix).dr = dr;
    ROIs(DMDix).fn = fn{DMDix};
    ROIs(DMDix).IM = IMt;
end
for DMDix = 1:nDMDs
    waitfor(hROIs(DMDix).hF);
    ROIs(DMDix).roiData = hROIs(DMDix).roiData;
end
save(fnAnn, 'ROIs', '-v7.3'); clear hROIs; 
disp('Post-Acquisition ROI annotations saved to:');
disp(fnAnn);


