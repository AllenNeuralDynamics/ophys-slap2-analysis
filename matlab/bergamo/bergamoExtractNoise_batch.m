dirs = uipickfiles();

dataDirs = {};
for i = 1:length(dirs)
    dataDirs = [dataDirs, findAllScanDirs(dirs{i})];
end

parfor idx = 1:length(dataDirs)
    [dFolder, dName] = fileparts(dataDirs{idx});
    disp(['Running directory ' fullfile(dFolder,dName)])

    try
        data = load(fullfile(dFolder,dName,[dName '_EXPTSUMMARY.mat']));
    catch
        try
            data = load(fullfile(dFolder,dName,[dName '_DENOISED_EXPTSUMMARY.mat']));
        catch
            disp(['Not processed: ' fullfile(dFolder,dName)]);
            continue;
        end
    end

    exptSummary = data.exptSummary;
    exptSummary.noiseEst = squeeze(estimatenoise(exptSummary.dFFls{1},2));

    for jdx = find(isnan(exptSummary.noiseEst))'
        [a,b] = ind2sub(size(exptSummary.noiseEst),jdx);
        try
            trace = exptSummary.dFFls{1}(a,:,b);
        catch
            disp(['Error estimating noise:' fullfile(dFolder,dName)])
            continue
        end
        times = 1:size(trace,2);
        if sum(isnan(trace)) < 0.9 * size(trace,2) % at least 10% of time points need to be valid to get a noise estimate
            times(isnan(trace)) = [];
            trace(isnan(trace)) = [];
            exptSummary.noiseEst(jdx) = estimatenoise(trace,times);
        end
    end

    savefunc(fullfile(dFolder,dName,[dName '_EXPTSUMMARY.mat']),exptSummary);
end

%%
function savefunc(fileName, exptSummary)
save(fileName,'exptSummary') %,'-v7.3') cannot read properly 7.3 in python using h5py
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