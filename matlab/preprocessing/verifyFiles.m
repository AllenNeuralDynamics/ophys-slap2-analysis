function [trialTable, keepTrials] = verifyFiles(fn,dr, params)

load([dr filesep fn], 'trialTable');
nTrials = length(trialTable.trueTrialIx);
nDMDs = size(trialTable.filename,1);
keepTrials = true(nDMDs, nTrials);
for trialIx = nTrials:-1:1
    for DMDix = 1:nDMDs
        [~, tiffFn] = fileparts(trialTable.fnRegDS{DMDix,trialIx}); 
        if isempty(dir([dr filesep tiffFn '.tif']))
            disp(['Missing tiff file:' tiffFn])
            keepTrials(DMDix,trialIx) = false;
        end
        [~, alignFn] = fileparts(trialTable.fnAdata{DMDix,trialIx}); 
        if ~exist([dr filesep alignFn '.mat'], 'file')
            disp(['Missing alignData file:' alignFn])
            keepTrials(DMDix,trialIx) = false;
        end
        sourceFn = trialTable.filename{DMDix,trialIx};
        if ~exist([dr filesep sourceFn], 'file')
            disp(['Missing source data file:' sourceFn])
            keepTrials(DMDix, trialIx) = false;
        end
        [~,~,ext] = fileparts(sourceFn);
        if strcmpi(ext, '.dat')
            trialTable.fnRaw{DMDix,trialIx} = sourceFn;
        else
            rawFn = trialTable.fnRaw{DMDix,trialIx};
            if ~exist([dr filesep rawFn], 'file')
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