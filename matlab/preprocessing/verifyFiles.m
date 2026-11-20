function [trialTable, keepTrials, dirs] = verifyFiles(fn, dr, params, dirs)
%VERIFYFILES Load a trial table and check that per-trial files exist.
%   [trialTable, keepTrials, dirs] = verifyFiles(fn, dr, params)
%   [trialTable, keepTrials, dirs] = verifyFiles(fn, dr, params, dirs)
%
%   fn may be trialTable.mat (legacy layout, everything in dr) or a
%   GIAnT-MATLAB trial_table.h5 (raw data in dirs.dataDr, registration
%   outputs in dirs.mocoDr). The table is returned in the legacy camelCase
%   layout via loadTrialTableLegacy. Filenames stay relative; dirs tells
%   callers where to resolve them.

if nargin < 4 || isempty(dirs)
    [trialTable, dirs] = loadTrialTableLegacy(fullfile(dr, fn));
else
    trialTable = loadTrialTableLegacy(fullfile(dr, fn));
end
nTrials = size(trialTable.filename, 2);
nDMDs = size(trialTable.filename,1);
if ~isfield(trialTable, 'fnRaw') || isempty(trialTable.fnRaw)
    trialTable.fnRaw = cell(nDMDs, nTrials);
end
keepTrials = true(nDMDs, nTrials);
for trialIx = nTrials:-1:1
    for DMDix = 1:nDMDs
        if isempty(trialTable.filename{DMDix,trialIx}) || isempty(trialTable.fnAdata{DMDix,trialIx})
            keepTrials(DMDix,trialIx) = false;
            continue
        end
        [~, tiffFn, ext] = fileparts(trialTable.fnRegDS{DMDix,trialIx});
        if isempty(dir(fullfile(dirs.mocoDr, [tiffFn '.tif']))) && isempty(dir(fullfile(dirs.mocoDr, [tiffFn '.h5'])))
            disp(['Missing tiff or h5 file:' tiffFn])
            keepTrials(DMDix,trialIx) = false;
        else
            trialTable.fnRegDS{DMDix,trialIx} = [tiffFn, ext];
        end

        [~, alignFn, alignExt] = fileparts(trialTable.fnAdata{DMDix,trialIx});
        if ~strcmpi(alignExt, '.h5')
            alignExt = '.mat';
        end
        if ~exist(fullfile(dirs.mocoDr, [alignFn alignExt]), 'file')
            disp(['Missing alignData file:' alignFn])
            keepTrials(DMDix,trialIx) = false;
        else
            trialTable.fnAdata{DMDix,trialIx} = [alignFn alignExt];
        end

        sourceFn = trialTable.filename{DMDix,trialIx};
        if ~exist(fullfile(dirs.dataDr, sourceFn), 'file')
            disp(['Missing source data file:' sourceFn])
            keepTrials(DMDix, trialIx) = false;
        end
        [~,~,ext] = fileparts(sourceFn);
        if strcmpi(ext, '.dat')
            trialTable.fnRaw{DMDix,trialIx} = sourceFn;
        else
            rawFn = trialTable.fnRaw{DMDix,trialIx};
            if isempty(rawFn) || ~exist(fullfile(dirs.dataDr, rawFn), 'file')
                disp(['Missing raw data file:' rawFn])
                keepTrials(DMDix, trialIx) = false;
            end
        end
    end
end
if ~all(keepTrials)
    disp(['Files were missing for ' int2str(sum(~keepTrials(:))) ' recordings; likely failed alignments. Proceeding without them.']);
end
if ~any(keepTrials)
    error('All trials were rejected due to missing alignment files!');
end
end
