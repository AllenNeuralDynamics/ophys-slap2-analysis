function trialTable = buildTrialTableBergamo(dr, fns)
%This function organizes multi-trial recordings and metadata for Bergamo
%recordings
%autoMode: no user input required

%get a list of tif files in a given folder
if ~nargin
    dr = uigetdir;
end
if nargin<2 
    unpickedfiles = dir([dr filesep '*.tif']);

    epoch = 0;
    while ~isempty(unpickedfiles)
        [indx,tf] = listdlg('ListString',{unpickedfiles.name}, 'PromptString',['Select files for EPOCH ' int2str(epoch)]);
        if ~tf
            break
        end
        epoch = epoch+1;
        epochfiles{epoch} = {unpickedfiles(indx).name};
        unpickedfiles(indx) = [];
    end
elseif ~iscell(fns) && fns==true %generate an autoTrialTable with all files
    %select all tif files in folder that are not REGISTERED and put them in
    %a single epoch
    epoch  = 1;
    files = dir([dr filesep '*.tif']);
    files = files(~contains({files.name}, 'REGISTERED'));
    epochfiles{1} = {files.name};
else %files were passed to generate an autoTrialTable
    epoch  = 1;
    epochfiles{1} = fns;
end

trialTable.filename = {};
trialTable.trialEndTimeFromPC = [];
trialTable.trialStartTimeInferred = [];

trueTrialIx = 0;
for eIx = 1:epoch %for each epoch
    files = epochfiles{eIx};
    for fIx = 1:length(files)
        trueTrialIx = trueTrialIx+1;
        trialTable.filename{1,trueTrialIx} = files{fIx};
        trialTable.trueTrialIx(trueTrialIx) = trueTrialIx;
        trialTable.epoch(trueTrialIx) = eIx;
    end
end

save([dr filesep 'trialTable'], 'trialTable');
end


