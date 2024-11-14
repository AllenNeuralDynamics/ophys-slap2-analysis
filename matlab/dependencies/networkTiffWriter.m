function networkTiffWriter(mat, fnwrite, pixelscale, localDr)

[dr, name, ext] = fileparts(fnwrite);

writeDr = dr;
isRemote = ~any(strcmpi(writeDr(1:2), {'C:', 'D:','E:', 'F:'}));
if isRemote
    if ~exist('localDr', 'var'); localDr = 'C:\temp'; end % 'F:\tmp_tiffIO'; end
  
    writeDr = localDr;
    if exist(writeDr,'dir') ~= 7; mkdir(writeDr); end
end

fTIF = Fast_BigTiff_Write(fullfile(writeDr, [name ext]),pixelscale,0);

for frame = 1:size(mat,3)
    fTIF.WriteIMG(single(mat(:,:,frame)));
end

fTIF.close;

if isRemote
    try
        copyfile(fullfile(writeDr, [name ext]),dr);
        delete(fullfile(writeDr, [name ext]));
    catch
        disp(['Could not copy file to dir ' dr])
    end
end

end