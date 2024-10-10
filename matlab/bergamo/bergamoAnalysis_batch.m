dirs = uipickfiles();

activityChannel = str2double(questdlg('Which color channel is glutamate activity in?', 'Set activityChannel','1','2','1'));
dataDirs = {};
for i = 1:length(dirs)
    dataDirs = [dataDirs, findAllScanDirs(dirs{i})];
end

parpool(8);
parfor idx = 1:length(dataDirs)
    [dFolder, dName] = fileparts(dataDirs{idx});
    disp(['Running directory ' dName])
    try
        % stripRegistrationBergamo_saveinplace([],fullfile(dFolder,dName,[dName '.tif']));
        %stripRegistrationBergamo([],[],fullfile(dFolder,dName,[dName '.tif']));
    
        % if ~exist(fullfile(dFolder,dName,[dName '_DENOISED_ALIGNMENTDATA.mat']),'file')
        %     copyfile(fullfile(dFolder,dName,[dName '_ALIGNMENTDATA.mat']),fullfile(dFolder,dName,[dName '_DENOISED_ALIGNMENTDATA.mat']));
        % end

        % tmp = dir(fullfile(dFolder,dName,[dName '_REGISTERED_RAW*']));
        % registeredRawFile = tmp.name;
        % denoise20Hz(fullfile(dFolder,dName),registeredRawFile);

        % tmp = dir(fullfile(dFolder,dName,'*_DENOISED_REGISTERED_DOWNSAMPLED-2x*'));
        tmp = dir(fullfile(dFolder,dName,'*_REGISTERED_DOWNSAMPLED-2x*'));
        downsampledFile = tmp.name;

        summarizeBergamo_Peaks(fullfile(dFolder,dName),downsampledFile, activityChannel);
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