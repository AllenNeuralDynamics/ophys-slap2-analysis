
params.doRegistration = false;
params.doDenoise = false;

sParams = setParams('summarizeBergamo_NoLoCo');
dirs = uipickfiles();

dataDirs = {};
for i = 1:length(dirs)
    dataDirs = [dataDirs, findScanDirsRecursive(dirs{i})];
end

%initialize parallel pool with n workers
%parpool(10);

parfor idx = 1:length(dataDirs)
    [dFolder, dName] = fileparts(dataDirs{idx});
    disp(['Running directory ' dName])
    try
        % if params.doRegistration
        %     stripRegistrationBergamo_saveinplace([],fullfile(dFolder,dName,[dName '.tif']));
        %     %stripRegistrationBergamo([],[],fullfile(dFolder,dName,[dName '.tif']));
        % end
        % 
        % if params.doDenoise
        %     % tmp = dir(fullfile(dFolder,dName,[dName '_REGISTERED_RAW*']));
        %     % registeredRawFile = tmp.name;
        %     % denoise20Hz(fullfile(dFolder,dName),registeredRawFile);
        %             % tmp = dir(fullfile(dFolder,dName,[dName '_REGISTERED_RAW*']));
        %     % registeredRawFile = tmp.name;
        %  % denoise20Hz(fullfile(dFolder,dName),registeredRawFile);
        % 
        %          % if ~exist(fullfile(dFolder,dName,[dName '_DENOISED_ALIGNMENTDATA.mat']),'file')
        % %     copyfile(fullfile(dFolder,dName,[dName '_ALIGNMENTDATA.mat']),fullfile(dFolder,dName,[dName '_DENOISED_ALIGNMENTDATA.mat']));
        % % end
        % end
        % 
        %copy locally if no file found
        % if ~exist(fullfile(dFolder,dName,[dName '_DENOISED_ALIGNMENTDATA.mat']),'file')
        %     disp('FLAG 1')
        %     copyfile(fullfile(dFolder,dName,[dName '_ALIGNMENTDATA.mat']),fullfile(dFolder,dName,[dName '_DENOISED_ALIGNMENTDATA.mat'])); %inputs are alignment and denoised alignment mat files
        % end

        %run processing
        tmp = dir([dataDirs{idx} filesep '*_DENOISED_REGISTERED_DOWNSAMPLED-2x*.tif']);
        if isempty(tmp)
            tmp = dir([dataDirs{idx} filesep '*_REGISTERED_DOWNSAMPLED-2x*.tif']);
        end
        downsampledFile = tmp(1).name;
        summarizeBergamo_NoLoCo(downsampledFile,dataDirs{idx}, sParams);
    catch ME
        fprintf("%s: %s\n",dName, ME.message)
    end
end

%%
function scanDirs = findAllScanDirs(dirName)

scanDirs = {};

[~, dirNameOnly] = fileparts(dirName);
if (length(dirNameOnly) > 5) && all(dirNameOnly(1:5) == 'scan_')
    scanDirs = [scanDirs, dirName];
    return;
end

subDirs = dir(dirName);
subDirs = subDirs([subDirs.isdir] & ~ismember({subDirs.name},{'.','..'}));

for idx = 1:length(subDirs)
    subDirPath = fullfile(dirName, subDirs(idx).name);
    scanDirs = [scanDirs, findAllScanDirs(subDirPath)];
end

end

%%
function scanTiffs = findAllScanTiffs(dirName)

scanDirs = {};

[~, dirNameOnly] = fileparts(dirName);
if (length(dirNameOnly) > 5) && all(dirNameOnly(1:5) == 'scan_')
    scanDirs = [scanDirs, dirName];
    return;
end

subDirs = dir(dirName);
subDirs = subDirs([subDirs.isdir] & ~ismember({subDirs.name},{'.','..'}));

for idx = 1:length(subDirs)
    subDirPath = fullfile(dirName, subDirs(idx).name);
    scanDirs = [scanDirs, findAllScanDirs(subDirPath)];
end

end

function scanDirs = findScanDirsRecursive(dirName)

scanDirs = {};

subDirs = dir(dirName);
subDirs = subDirs([subDirs.isdir] & ~ismember({subDirs.name},{'.','..'}));

for idx = 1:length(subDirs)
    subDirPath = fullfile(dirName, subDirs(idx).name);
    scanDirs = [scanDirs, findScanDirsRecursive(subDirPath)];
end

dd = dir([dirName filesep '*_REGISTERED_DOWNSAMPLED-2x*.tif']);
containsTifs = ~isempty(dd);
if containsTifs
     scanDirs = [scanDirs dirName];
     return
end
end



