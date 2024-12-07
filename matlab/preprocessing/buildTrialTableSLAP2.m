function trialTable = buildTrialTableSLAP2(dr)
%This function organizes multi-trial recordings and reference images as a first step in the SLAP2 data processing
%pipeline

%This function addresses a bug in SLAP2 trial numbering as of Feb 2024, where trial
%numbers sometimes fail to increment. THis makes some files extra long,
% and subsequent trial numbers get out of sync

%parameters
lineDiffThresh = 2000; %difference threshold for calling two recordings the same length, in lines. ~0.2 seconds

%get a list of dat files in a given folder
if ~nargin
    dr = uigetdir;
end
unpickedfiles = dir([dr filesep '*.dat']);

epoch = 0;
while ~isempty(unpickedfiles)
    [indx,tf] = listdlg('ListString',{unpickedfiles.name}, 'PromptString',['Select files for EPOCH ' int2str(epoch)]);
    if ~tf
        break
    end
    epoch = epoch+1;
    epochfiles{epoch} = unpickedfiles(indx);
    unpickedfiles(indx) = [];
end

%first, load the reference images, deduce the imaging plane of the ROI, and
% %compute the soma ROI
DMDixs = [1 2];
%find the files within the entire folder structure named REFERENCE
for DMDix = DMDixs
    list = dir([dr filesep '**' filesep '*DMD' int2str(DMDix) '*-REFERENCE*']);
    switch length(list)
    case 0
        error(['Could not find reference image within folder: ' dr])
    case 1
        DMD1refFn = fullfile(list(1).folder, list(1).name);
        A = ScanImageTiffReader(DMD1refFn);
        refStack{DMDix}.IM = A.data;
    
        %load metadata for the reference stack and BCI ROI
        desc = A.descriptions;
        zIx = 0; Zs = []; channels = [];
        for imIx = 1:numel(desc)
                jj = jsondecode(desc{imIx});
                if ~any(channels == jj.channel)
                    channels = [channels jj.channel];
                end
                if imIx==1
                    accumChan = jj.channel;
                end
                if jj.channel == accumChan
                    zIx = zIx +1;
                    Zs(zIx) = jj.z;
                end
        end
        refStack{DMDix}.channels = channels;
        refStack{DMDix}.Zs = Zs;
        refStack{DMDix}.dmdPixel2SampleTransform = jj.dmdPixel2SampleTransform;

    % case 2  %old format, two reference files
    %     error('Too many reference stacks found in the specified directory!');
    %     % for cix = 1:2
    %     %     DMD1refFn = fullfile(list(cix).folder, list(cix).name);
    %     %     A = ScanImageTiffReader(DMD1refFn);
    %     %     refStack{DMDix}.IM(:,:,:,cix) = A.data;
    %     % end
    otherwise
        list = dir([dr filesep '**' filesep '*DMD' int2str(DMDix) '_CONFIG1-REFERENCE*']);
        switch length(list)
            case 0
                error('Too many reference stacks found in the specified directory!');
            case 1
                DMD1refFn = fullfile(list(1).folder, list(1).name);
                A = ScanImageTiffReader(DMD1refFn);
                refStack{DMDix}.IM = A.data;
            
                %load metadata for the reference stack and BCI ROI
                desc = A.descriptions;
                zIx = 0; Zs = []; channels = [];
                for imIx = 1:numel(desc)
                        jj = jsondecode(desc{imIx});
                        if ~any(channels == jj.channel)
                            channels = [channels jj.channel];
                        end
                        if imIx==1
                            accumChan = jj.channel;
                        end
                        if jj.channel == accumChan
                            zIx = zIx +1;
                            Zs(zIx) = jj.z;
                        end
                end
                refStack{DMDix}.channels = channels;
                refStack{DMDix}.Zs = Zs;
                refStack{DMDix}.dmdPixel2SampleTransform = jj.dmdPixel2SampleTransform;
            otherwise
                error('Too many CONFIG1 reference stacks found in the specified directory!');
        end
    end


end
trialTable.refStack = refStack; clear refStack;

trialTable.filename = {};
trialTable.firstLine = [];
trialTable.lastLine = [];
trialTable.trialEndTimeFromPC = [];
trialTable.trialStartTimeInferred = [];

trueTrialIx = 0;
for eIx = 1:epoch %for each epoch
    %get number of lines from dat file timestamp
    files = epochfiles{eIx};
    disp(['Loading metadata from ' int2str(length(epochfiles{eIx})) ' DAT files...']);
    for fIx = length(files):-1:1
        hDat = slap2.Slap2DataFile([dr filesep files(fIx).name]);
        numLines(fIx) = hDat.totalNumLines;
    end
    linePeriod_s = hDat.hDataFile.metaData.linePeriod_s;
    disp('done loading ')
    isDMD1 = cellfun(@(x)(~isempty(x)), strfind({files.name}, 'DMD1'));

    DMD1files = files(isDMD1);
    [~, sortorder] = sort([DMD1files.datenum], 'ascend');
    DMD1files = DMD1files(sortorder);
    numLinesDMD1 = numLines(isDMD1); numLinesDMD1 = numLinesDMD1(sortorder);

    DMD2files = files(~isDMD1);
    [~, sortorder] = sort([DMD2files.datenum], 'ascend');
    DMD2files = DMD2files(sortorder);
    numLinesDMD2 = numLines(~isDMD1); numLinesDMD2 = numLinesDMD2(sortorder);

    lastDMD1fIx = 0;
    lastDMD2fIx = 0;
    accumLines = [0 0];
    while lastDMD2fIx<length(DMD2files) && lastDMD1fIx<length(DMD1files)

        %ocnfirm that the trials match up AT SOME POINT SOON
        cumLines1 = cumsum(numLinesDMD1(lastDMD1fIx+1: min(end, lastDMD1fIx+3)))-accumLines(1);
        cumLines2 = cumsum(numLinesDMD2(lastDMD2fIx+1: min(end, lastDMD2fIx+3)))-accumLines(2);
        assert(min(abs(cumLines1 - cumLines2'), [],'all')<lineDiffThresh, 'Error lining up trials across DMDs!')

        nLines = min(cumLines1(1), cumLines2(1));
        trueTrialIx = trueTrialIx+1;
        trialTable.filename{1,trueTrialIx} = DMD1files(lastDMD1fIx+1).name;
        trialTable.filename{2,trueTrialIx} = DMD2files(lastDMD2fIx+1).name;
        trialTable.firstLine(1,trueTrialIx) = accumLines(1)+1;
        trialTable.firstLine(2,trueTrialIx) = accumLines(2)+1;
        trialTable.lastLine(1,trueTrialIx) = accumLines(1)+nLines;
        trialTable.lastLine(2,trueTrialIx) = accumLines(2)+nLines;
        trialTable.trueTrialIx(trueTrialIx) = trueTrialIx;
        trialTable.epoch(trueTrialIx) = eIx;

        if abs(cumLines1(1)-cumLines2(1))<lineDiffThresh
            trialTable.trialEndTimeFromPC(trueTrialIx) = DMD1files(lastDMD1fIx+1).datenum;
            trialTable.trialStartTimeInferred(trueTrialIx) = trialTable.trialEndTimeFromPC(trueTrialIx) - datenum(duration(0,0,nLines.*linePeriod_s));

            lastDMD1fIx = lastDMD1fIx+1; %we finished this DMD1 file
            lastDMD2fIx = lastDMD2fIx+1; %we finished this DMD2 file
            accumLines = [0 0]; %reset accumulated lines

        elseif cumLines1(1) < cumLines2(1)
            trialTable.trialEndTimeFromPC(trueTrialIx) = DMD1files(lastDMD1fIx+1).datenum;
            trialTable.trialStartTimeInferred(trueTrialIx) = trialTable.trialEndTimeFromPC(trueTrialIx) - datenum(duration(0,0,nLines.*linePeriod_s));

            lastDMD1fIx = lastDMD1fIx+1; %we finished this DMD1 file
            accumLines(1) = 0;
            accumLines(2) = accumLines(2) + nLines;
        elseif cumLines2(1) < cumLines1(1)
            trialTable.trialEndTimeFromPC(trueTrialIx) = DMD2files(lastDMD2fIx+1).datenum;
            trialTable.trialStartTimeInferred(trueTrialIx) = trialTable.trialEndTimeFromPC(trueTrialIx) - datenum(duration(0,0,nLines.*linePeriod_s));

            lastDMD2fIx = lastDMD2fIx+1; %we finished this DMD1 file
            accumLines(1) = accumLines(1) + nLines;
            accumLines(2) = 0;
        end
    end
end


save([dr filesep 'trialTable'], 'trialTable');
end


