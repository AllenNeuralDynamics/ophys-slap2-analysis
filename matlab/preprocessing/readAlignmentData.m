function aData = readAlignmentData(path)
%READALIGNMENTDATA Load *_ALIGNMENTDATA.mat (multiRoiRegSLAP2) or
%   *_ALIGNMENTDATA.h5 (GIAnT-MATLAB MultiRoiRegistration) as the legacy
%   flat aData struct expected by loadAndProcessVoltage.
%
%   GIAnT nests SLAP2-specific fields under aData.slap2 with different
%   names; they are hoisted here:
%     slap2.onlineMotionXshift/Yshift/Zshift -> onlineXshift/Yshift/Zshift
%     slap2.Z_depths                         -> Z
%     slap2.varFacDS, cropRow, cropCol, viewC, viewR, trimRows, trimCols -> top level
%   aError is not saved by GIAnT; it is filled with NaN so downstream
%   upsampling still works.

[~, ~, ext] = fileparts(char(path));
if strcmpi(ext, '.h5')
    aData = readStructH5(path);
    if isfield(aData, 'slap2') && isstruct(aData.slap2)
        sl = aData.slap2;
        renames = {'onlineMotionXshift', 'onlineXshift'; ...
                   'onlineMotionYshift', 'onlineYshift'; ...
                   'onlineMotionZshift', 'onlineZshift'; ...
                   'Z_depths', 'Z'};
        slFields = fieldnames(sl);
        for ix = 1:numel(slFields)
            name = slFields{ix};
            hit = strcmp(renames(:,1), name);
            if any(hit)
                name = renames{hit, 2};
            end
            if ~isfield(aData, name)
                aData.(name) = sl.(slFields{ix});
            end
        end
        aData = rmfield(aData, 'slap2');
    end
    if ~isfield(aData, 'aError') && isfield(aData, 'recNegErr')
        aData.aError = nan(size(aData.recNegErr));
    end
    vecFields = {'onlineXshift', 'onlineYshift', 'onlineZshift', 'motionDSc', 'motionDSr', 'recNegErr', 'aError'};
    for ix = 1:numel(vecFields)
        if isfield(aData, vecFields{ix})
            aData.(vecFields{ix}) = double(aData.(vecFields{ix})(:).');
        end
    end
else
    S = load(path, 'aData');
    aData = S.aData;
end
end
