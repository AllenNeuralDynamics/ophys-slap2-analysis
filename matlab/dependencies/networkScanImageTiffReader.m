function [data, desc, meta] = networkScanImageTiffReader(fname, localDr)

import ScanImageTiffReader.ScanImageTiffReader;

[dr, name, ext] = fileparts(fname);
tmpFileName = [name ext];
isRemote = ~any(strcmpi(dr(1:2), {'C:', 'D:','E:', 'F:'}));

if isRemote
    if ~exist('localDr', 'var'); localDr = 'C:\temp'; end % 'F:\tmp_tiffIO'; end

    dr = localDr;

    if exist(dr,'dir') ~= 7; mkdir(dr); end

    try
        tmpFileName = [name '_' sprintf('%03d',randi(999)) ext];
        copyfile(fname,fullfile(dr,tmpFileName));
    catch
        disp(['Could not copy file ' fname]);
        data = nan;
        desc = nan;
        meta = nan;
    end
end


A = ScanImageTiffReader(fullfile(dr,tmpFileName));
desc = A.descriptions();
meta = A.metadata;

data = single(A.data);

clear('A');

if isRemote
    delete(fullfile(localDr, tmpFileName));
end

end