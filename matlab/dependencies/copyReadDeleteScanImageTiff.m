function [data, meta]= copyReadDeleteScanImageTiff(remotepath, localDir)
if isempty(remotepath) %PASS EMPTY MATRIX TO INITIALIZE
    thisDir = fileparts(which('copyReadDeleteScanImageTiff'));
    if ~exist([thisDir filesep 'pathToFastDrive.mat'], 'file')
        localDir = uigetdir('', 'Select a local directory to save temporary files, on your FASTEST local drive!');
        if ~exist(localDir, 'dir')
            warning('User canceled setting local directory. defaulting to C:/temp')
            localDir = 'C:\temp\';
        end
        save([thisDir filesep 'pathToFastDrive.mat'], 'localDir');
    end
    return
end

if nargin<2
    %read user settings for local directory
    thisDir = fileparts(which('copyReadDeleteScanImageTiff'));
    if exist([thisDir filesep 'pathToFastDrive.mat'], 'file')
        load([thisDir filesep 'pathToFastDrive.mat'], 'localDir');
    else
        localDir = uigetdir('', 'Select a local directory to save temporary files, on your FASTEST local drive!');
        if exist(localDir, 'dir')
           save([thisDir filesep 'pathToFastDrive.mat'], 'localDir');
        else
            warning('User canceled setting local directory. defaulting to C:/temp')
            localDir = 'C:\temp\';
        end
    end
end
if ~exist(localDir, 'dir')
    mkdir(localDir);
end

assert(exist(remotepath, 'file'));
randName = [int2str(round(1e10*rand+0.1)) '.tif'];
localpath = [localDir  filesep   randName];
disp(localpath)
copyfile(remotepath, localpath);
A = ScanImageTiffReader(localpath);
data = A.data();
if nargout>1
    meta = A.metadata;
end
delete(A);
delete(localpath);
end