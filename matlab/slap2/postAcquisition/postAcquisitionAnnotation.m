function postAcquisitionAnnotation (dr)

if ~nargin
    dr = pwd;
end
import ScanImageTiffReader.ScanImageTiffReader;

fn{1} = uigetfile('*DMD1*.tif', 'Select a DMD1 acquisition file for ANNOTATION');
fn{2} = uigetfile('*DMD2*.tif', 'Select a DMD2 acquisition file for ANNOTATION');

fnAnn = [dr filesep 'ANNOTATIONS.mat'];

if exist(fnAnn, 'file')
    load(fnAnn, 'ROIs')
else
    for DMDix = 1:nDMDs
        %load image data
        A = ScanImageTiffReader(fullfile(dr,fn{DMDix}));
        IM  = A.data;
        desc = A.descriptions();
        meta = jsondecode(desc{1});
        nchan = meta.numChannels;
        nframes = size(IM,3)./nchan;

        IM = double(reshape(IM, size(IM,1), size(IM,2), nchan, nframes));
        IMt = squeeze(trimmean(IM, 40, 4));

        hROIs(DMDix) = drawROIs(sqrt(max(0,IMt(:,:,end:-1:1))), dr, fn{DMDix}); %3rd dim inverted because channels are GR and colors are RGB
        
        ROIs(DMDix).dr = dr;
        ROIs(DMDix).fn = fn{DMDix};
        ROIS(DMDix.IM) = IM
    end
    for DMDix = 1:nDMDs
        waitfor(hROIs(DMDix).hF);
        ROIs(DMDix).roiData = hROIs(DMDix).roiData;
    end
    save(fnAnn, 'ROIs'); clear hROIs;
end

