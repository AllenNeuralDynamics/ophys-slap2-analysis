function s = readStructH5(filename)
%READSTRUCTH5 Read an HDF5 file written by GIAnT-MATLAB saveStructToH5.
%   Groups become nested structs; string datasets become cellstr (filename
%   grids keep their shape), scalar strings become char. Mirrors
%   GIAnT-MATLAB/dependencies/io/loadStructFromH5 without requiring it on
%   the path.
if ~exist(filename, 'file')
    error('readStructH5:FileNotFound', 'File not found: %s', filename);
end
info = h5info(filename);
s = readGroup(filename, info);
end

function s = readGroup(filename, grpInfo)
s = struct();
cellFields = {'filename', 'fn_reg_ds', 'fn_raw', 'fn_adata', 'fn_ann'};
for ix = 1:numel(grpInfo.Datasets)
    ds = grpInfo.Datasets(ix);
    if strcmp(grpInfo.Name, '/')
        path = ['/' ds.Name];
    else
        path = [grpInfo.Name '/' ds.Name];
    end
    val = h5read(filename, path);
    if isstring(val)
        if any(strcmp(ds.Name, cellFields))
            val = reshape(cellstr(val(:)), size(val));
        elseif isscalar(val)
            val = char(val);
        else
            val = cellstr(val);
        end
    elseif ischar(val) && any(strcmp(ds.Name, cellFields)) && isvector(val)
        val = {strtrim(val(:).')};
    end
    s.(ds.Name) = val;
end
for ix = 1:numel(grpInfo.Groups)
    grp = grpInfo.Groups(ix);
    parts = strsplit(grp.Name, '/');
    s.(parts{end}) = readGroup(filename, grp);
end
end
