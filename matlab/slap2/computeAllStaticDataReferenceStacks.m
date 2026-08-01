function computeAllStaticDataReferenceStacks(rootDir)

if nargin < 1 || isempty(rootDir)
    rootDir = uigetdir(pwd, 'Select root directory');
    if isequal(rootDir, 0)
        return;
    end
end

allItems = dir(fullfile(rootDir, '**'));

staticDirs = allItems([allItems.isdir] & strcmp({allItems.name}, 'static_data'));
staticDirPaths = fullfile({staticDirs.folder}, {staticDirs.name});

missingReference = {};

for d = 1:numel(staticDirPaths)

    folder = staticDirPaths{d};
    tifFiles = dir(fullfile(folder, '*.tif'));

    for k = 1:numel(tifFiles)

        name = tifFiles(k).name;

        if endsWith(name, '-REFERENCE.tif', 'IgnoreCase', true)
            continue
        end

        [~, base, ext] = fileparts(name);
        refFile = fullfile(folder, [base '-REFERENCE' ext]);

        if ~isfile(refFile)
            missingReference{end+1,1} = fullfile(folder, name);
        end
    end
end

for ix = 1:numel(missingReference)
    disp(missingReference{ix});
    obj = slap2.gui.refstack.ReferenceStack.loadTif(missingReference{ix});
    clear obj;
end