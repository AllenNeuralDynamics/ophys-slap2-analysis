function stripRegBergamo(dr, fns, paramsIn)
if ~nargin || isempty(dr)
    [fns, dr] = uigetfile('*.*', 'Select either a trialTable.mat file or your tifs to register', 'MultiSelect','on');
end
if iscell(fns) %user selected multiple tiffs; generate a trial table
    trialTable = buildTrialTableBergamo(dr, fns);
elseif contains(fns, '.tif') %user selected single tiff; generate a trial table
    trialTable = buildTrialTableBergamo(dr, fns);
elseif contains(fns, 'trialTable')
    load([dr filesep fns], 'trialTable');
else
    error('Must select either TIF files or a trial table file')
end
fullPathToTrialTable = [dr filesep 'trialTable.mat'];

%PARAMETER SETTING
if nargin>2
    params = setParams('stripRegBergamo', paramsIn);
else
    params = setParams('stripRegBergamo');
end

%set up parallelization
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
nWorkers = min(params.nWorkers, numel(trialTable.filename));
if poolsize<nWorkers
    delete(gcp('nocreate'));
    if params.nWorkers<4
        warning('You are using few parallel workers!');
    end
    parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM %4-30-24, lowering processes again to prevent another error (18 --> 15)
end

%align in parallel
nTrials= numel(trialTable.filename);
fnRegDS = cell(1, nTrials); fnRaw = cell(1, nTrials); fnAdata = cell(1, nTrials);
parfor f_ix = 1:nTrials
    [fnRegDS{f_ix}, fnRaw{f_ix}, fnAdata{f_ix}]= alignAsync(dr, trialTable, params, f_ix);
end

trialTable.fnRegDS = fnRegDS;
trialTable.fnRaw = fnRaw;
trialTable.fnAdata = fnAdata;
trialTable.alignParams = params;
save(fullPathToTrialTable , "trialTable")

disp('done bergamoRegistration.')

end


function [fnDS, fnRaw, fnAdata] = alignAsync(dr, trialTable, params, f_ix);
maxshift = params.maxshift;
dsFac = 2^params.ds_time; params.dsFac = dsFac;

fn = trialTable.filename{f_ix};
[~,fn, ext] = fileparts(fn);
fnstem = fn;
fn = strcat(fn,ext);

disp(['Aligning: ' [dr filesep fn]])

[Ad, desc, meta] = networkScanImageTiffReader([dr filesep fn]);
evalc(desc{1});

metaLines = strsplit(meta, '\n');
for lineIx = 1:length(metaLines)
    try
        evalc([metaLines{lineIx} ';']);
    catch
    end
end
numChannels = length(SI.hChannels.channelSave);
selCh = 1:numChannels;

pat = "frameTimestamps_sec = " + digitsPattern + "." + digitsPattern;
for frame = 1:10*numChannels %compute the framerate from the metadata by reading a few frames
    E = extract(desc{frame}, pat);
    timestamp(frame) = str2double(E{1}(23:end)); %#ok<AGROW>
end

frametime = median(diff(timestamp(1:numChannels:end)));
aData.numChannels = numChannels;
aData.frametime = frametime;
aData.alignHz = 1/frametime/dsFac;

Ad = permute(reshape(Ad, size(Ad,1), size(Ad,2), numChannels, []), [2 1 3 4]);
Ad = Ad(params.removeLines+1:end,:,:,:);

%subtract baseline
tsamp = round(linspace(1,size(Ad,4), 500));
BG = min(1000,prctile(Ad(:,:,:,tsamp), 10, [1 2 4]));
Ad = Ad-min(BG,1000);

%make an initial template with normcorr
initFrames = 1000;
framesToRead = initFrames * dsFac;
try
    Y = downsampleTime(Ad(:,:,:,1:framesToRead), params.ds_time);
catch ME
    disp(['Your file was too short:' fn])
    fnDS = []; fnRaw=[]; fnAdata = [];
    return
end
sz = size(Ad);
Yhp = squeeze(sum(reshape(Y(:,:,selCh,:),size(Y,1),size(Y,2),[],size(Y,4)),3));

% find the most correlated ~100 frames within the first initFrames to
% use as an initial template
rho = corr(reshape(Yhp,[],size(Yhp,3)));
dist_matrix = 1 - rho;
Z = linkage(squareform(dist_matrix),'average');

cutoff = 0.01;
min_cluster_size = 100;
clusters = [];
max_cutoff = 2.0;
while isempty(clusters) || all(cellfun(@numel, clusters) < min_cluster_size)
    cutoff = cutoff + 0.01;
    if cutoff > max_cutoff
        error('Cound not find a cluster with at least %d samples',min_cluster_size);
    end
    T = cluster(Z, 'cutoff', cutoff, 'criterion', 'distance');
    clusters = arrayfun(@(x) find(T == x), unique(T), 'UniformOutput', false);
end

max_mean_corr = -Inf;
for i = 1:length(clusters)
    cluster_indices = clusters{i};
    if numel(cluster_indices) >= min_cluster_size
        mean_corr = mean(mean(rho(cluster_indices,cluster_indices)));
        if mean_corr > max_mean_corr
            max_mean_corr = mean_corr;
            best_cluster = cluster_indices;
        end
    end
end

options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',initFrames,'max_shift',maxshift,'us_fac',4,'init_batch',initFrames, 'correct_bidir', false);
F = normcorre(Yhp,options_rigid,mean(Yhp(:,:,best_cluster),3));
F = mean(F,3); %fixed image
% F = mean(sqrt(abs(F)).*sign(F),3); %fixed image
template = nan(2*maxshift+sz(1), 2*maxshift+sz(2));
templateFull = template;
T0 = template; T00 = zeros(size(template)); templateCt = zeros(size(template));
T0(maxshift+(1:sz(1)), maxshift+(1:sz(2)))=F;
clear Y Yhp;

%aligned = nan(2*maxshift+szY(1), 2*maxshift+szY(2), nDSframes);
%alignedY = nan(2*maxshift+szY(1), 2*maxshift+szY(2), nDSframes);
initR = 0; initC = 0;
nDSframes= floor(sz(4)./(dsFac)); %number of downsampled frames
motionDSr = nan(1,nDSframes); motionDSc = nan(1,nDSframes); %matrices to store the inferred motion
aErrorDS = nan(1,nDSframes);
[viewR, viewC] = ndgrid((1:(sz(1)+2*maxshift))-maxshift, (1:(sz(2)+2*maxshift))-maxshift); %view matrices for interpolation

%for spatial downsampling, used to calculate alignment quality and for quicker visualization
dsTimes = 2;
dsSz = floor(size(template)./(2^dsTimes));
A_ds = nan([dsSz nDSframes]);

disp('Registering:');
for DSframe = 1:nDSframes
    readFrames = (DSframe-1)*(dsFac) + (1:(dsFac));

    M = downsampleTime(Ad(:,:,:, readFrames), params.ds_time);
    M = squeeze(sum(reshape(M(:,:,selCh,:),size(M,1),size(M,2),[],size(M,4)),3)); %merge colors

    if ~mod(DSframe, 1000)
        disp([int2str(DSframe) ' of ' int2str(nDSframes)]);
    end

    Ttmp = mean(cat(3, T0, T00,template),3, 'omitnan');
    % T = Ttmp(maxshift-initR + (1:sz(1)), maxshift-initC+(1:sz(2)));
    T = Ttmp(maxshift-initR + (1:sz(1)), maxshift-initC+(1:sz(2)));

    % should edit template so that it doesn't just use the cropped area
    % want to use the additional information in the template around the
    % init crop area

    % Mfull = interp2(1:sz(2), 1:sz(1), M,viewC, viewR, 'linear', nan);

    %output = dftregistration(fft2(M),fft2(single(T)),4);
    output = dftregistration_clipped(fft2(M),fft2(single(T)),4, params.clipShift);
    % [motion, R] = xcorr2_nans(Mfull,Ttmp,[initR;initC],params.clipShift);

    %         if abs(output(3))>10 || abs(output(4))>10
    %             %the shift was very large- what's up?
    %         end
    motionDSr(DSframe) = initR+output(3);
    motionDSc(DSframe) = initC+output(4);
    aErrorDS(DSframe) = output(1);

    % motionDSr(DSframe) = motion(1);
    % motionDSc(DSframe) = motion(2);
    % aErrorDS(DSframe) = R;

    if sqrt((motionDSr(DSframe)/sz(1)).^2 + (motionDSc(DSframe)/sz(2)).^2) > 0.75^2
        Mfull = interp2(1:sz(2), 1:sz(1), M,viewC, viewR, 'linear', nan);
        [motion, R] = xcorr2_nans(Mfull,Ttmp,[initR;initC],50);

        motionDSr(DSframe) = motion(1);
        motionDSc(DSframe) = motion(2);
        aErrorDS(DSframe) = R;
    end

    if abs(motionDSr(DSframe))<maxshift && abs(motionDSc(DSframe))<maxshift
        A = interp2(1:sz(2), 1:sz(1), M,viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan);
        % sel = ~isnan(A);

        %downsample in space and save, we will compute error statistics
        %from the downsampled data
        dsTmp = A;
        for dsIx = 1:dsTimes
            dsTmp = dsTmp(1:2:2*floor(end/2), 1:2:2*floor(end/2)) + dsTmp(1:2:2*floor(end/2), 2:2:2*floor(end/2)) + dsTmp(2:2:2*floor(end/2), 1:2:2*floor(end/2)) + dsTmp(2:2:2*floor(end/2), 2:2:2*floor(end/2));
        end
        A_ds(:,:,DSframe) = dsTmp;

        % Asmooth = imgaussfilt(A,1);
        %
        % selCorr = ~(isnan(Asmooth) | isnan(Ttmp));
        % aRankCorr(DSframe) = corr(Asmooth(selCorr), Ttmp(selCorr), 'type', 'Spearman');
        % recNegErr(DSframe) = mean(min(0, Asmooth(selCorr) .* mean(Ttmp(selCorr)) ./ mean(Asmooth(selCorr)) - Ttmp(selCorr)).^2);

        templateFull = sum(cat(3,templateFull .* templateCt, A),3,"omitnan");
        templateCt = templateCt + ~isnan(A);
        templateFull = templateFull ./ templateCt;
        template = templateFull;
        template(templateCt < 100) = nan;

        initR = round(motionDSr(DSframe));
        initC = round(motionDSc(DSframe));
    else
        motionDSr(DSframe) = initR;
        motionDSc(DSframe) = initC;
    end
end

%compute alignemnt error stats
offset = max(0, prctile(template(:), 25));
recNegErr = nan(1,nDSframes);
nChunks = ceil(nDSframes./(aData.alignHz*10)); %align 10 seconds at a time
chunkEdges = round(linspace(1, nDSframes+1, nChunks));
for chunkIx = 1:length(chunkEdges)-1
    t_ixs = chunkEdges(chunkIx):chunkEdges(chunkIx+1)-1;
    template = median(A_ds(:,:,t_ixs), 3, 'omitmissing');
    nanFrac = mean(isnan(A_ds(:,:,t_ixs)), 3);
    template(nanFrac>0.2) = nan;
    template_gamma = sqrt(max(0,template)+offset);
    template = repmat(template, 1,1,length(t_ixs));
    template(isnan(A_ds(:,:,t_ixs))) = nan;

    recNegErr(1,t_ixs) = sqrt(squeeze(mean((max(0, (template-A_ds(:,:,t_ixs))./template_gamma).^2), [1 2], 'omitnan')./mean(max(0,(template./template_gamma).^2, 'includenan'), [1 2], 'omitnan')));
end

%upsample the shifts and compute a tighter field of view
tDS = (1:nDSframes).*(dsFac) - (2^(params.ds_time-1)) + 0.5;
motionC = interp1(tDS, motionDSc, 1:((2^params.ds_time)*nDSframes), 'pchip','extrap'); %upsample the movement to the original timebase
motionR = interp1(tDS, motionDSr, 1:((2^params.ds_time)*nDSframes), 'pchip', 'extrap'); %upsample the movement to the original timebase
%aError = interp1(tDS, aErrorDS, 1:((2^params.ds_time)*nDSframes), 'nearest', 'extrap');
maxshiftC = max(abs(motionC)); maxshiftR = max(abs(motionR));
[viewR, viewC] = ndgrid((1:(sz(1)+2*maxshiftR))-maxshiftR, (1:(sz(2)+2*maxshiftC))-maxshiftC); %view matrices for interpolation

pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM

%save a downsampled aligned recording
fnDS = [fnstem '_REGISTERED_DOWNSAMPLED-' int2str(dsFac) 'x.tif'];
% fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
Bsum = zeros([size(viewR') numChannels]);
Bcount = zeros([size(viewR') numChannels]);

tiffSave = single(zeros([size(viewR, [2 1]) nDSframes*numChannels]));

for DSframe = 1:nDSframes
    readFrames = (DSframe-1)*(dsFac) + (1:(dsFac));
    YY = downsampleTime(Ad(:,:,:, readFrames), params.ds_time);
    for ch = 1:numChannels
        B =  interp2(1:sz(2), 1:sz(1), YY(:,:,ch),viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan)';
        tiffSave(:,:,(DSframe-1)*numChannels+ch) = single(B);
        % fTIF.WriteIMG(single(B));

        Bcount(:,:,ch) = Bcount(:,:,ch) + ~isnan(B);
        B(isnan(B)) = 0;
        Bsum(:,:,ch) = Bsum(:,:,ch)+double(B);
    end
end

nanRows = mean(mean(isnan(tiffSave),3),2) == 1;
nanCols = mean(mean(isnan(tiffSave),3),1) == 1;

tiffSave(nanRows,:,:) = [];
tiffSave(:,nanCols,:) = [];

networkTiffWriter(tiffSave, [dr filesep fnDS], pixelscale);
clear('tiffSave');

% %save an average image for each channel
% for ch = 1:numChannels
%     Bmean = Bsum(:,:,ch)./Bcount(:,:,ch);
%     minV = prctile(Bmean(~isnan(Bmean(:))), 10);
%     maxV = prctile(Bmean(~isnan(Bmean(:))), 99.9);
%     Bmean = uint8(255*sqrt(max(0,(Bmean-minV)./(maxV-minV))));
%     fnwrite = [fnstem '_REGISTERED_AVG_CH' num2str(ch) '_8bit.tif'];
%     fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
%     fTIF.WriteIMG(single(Bmean));
%     fTIF.close;
% end

%save an average image for each channel
for ch = 1:numChannels
    Bmean = Bsum(:,:,ch)./Bcount(:,:,ch);
    minV = prctile(Bmean(~isnan(Bmean(:))), 10);
    maxV = prctile(Bmean(~isnan(Bmean(:))), 99.9);
    Bmean = uint8(255*sqrt(max(0,(Bmean-minV)./(maxV-minV))));
    Bmean(nanRows,:) = [];
    Bmean(:,nanCols) = [];
    fnwrite = [fnstem '_REGISTERED_AVG_CH' num2str(ch) '_8bit.tif'];
    networkTiffWriter(single(Bmean), [dr filesep fnwrite], pixelscale);
    % fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
    % fTIF.WriteIMG(single(Bmean));
    % fTIF.close;
end

% %save an original-time-resolution recording
% fnwrite = [fnstem '_REGISTERED_RAW.tif'];
% fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
% for frame = 1:length(motionC)
%     for ch = 1:numChannels
%         B = interp2(1:sz(2), 1:sz(1), Ad(:,:,ch,frame),viewC+motionC(frame), viewR+motionR(frame), 'linear', nan)';
%         fTIF.WriteIMG(single(B));
%     end
% end
% fTIF.close;

fnRaw = [fnstem '_REGISTERED_RAW.tif'];
% fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
tiffSave = single(zeros([size(viewR, [2 1]) - [sum(nanRows) sum(nanCols)] length(motionC)*numChannels]));
for frame = 1:length(motionC)
    for ch = 1:numChannels
        B = interp2(1:sz(2), 1:sz(1), Ad(:,:,ch,frame),viewC+motionC(frame), viewR+motionR(frame), 'linear', nan)';
        B(nanRows,:) = [];
        B(:,nanCols) = [];
        tiffSave(:,:,(frame-1)*numChannels+ch) = single(B);
        % fTIF.WriteIMG(single(B));
    end
end

networkTiffWriter(single(tiffSave), [dr filesep fnRaw], pixelscale);
clear('tiffSave')

%save alignment metadata
aData.motionC =  motionC;
aData.motionR = motionR;
aData.motionDSc = motionDSc;
aData.motionDSr = motionDSr;
%aData.aError = aError;
%aData.aRankCorr = aRankCorr;
aData.recNegErr = recNegErr;
fnAdata = [fnstem '_ALIGNMENTDATA.mat'];
save([dr filesep fnAdata], 'aData');
end


function Y = downsampleTime(Y, ds_time)
for ix = 1:ds_time
    Y = Y(:,:,:,1:2:(2*floor(end/2)))+ Y(:,:,:,2:2:end);
end
Y = Y./2.^ds_time;
end

function [Y, readsuccess, done] = readFrames(tiffObj, nFramesToRead, nChannels)
nReads = nFramesToRead*nChannels;
done = false;
Y = single(tiffObj.read);
Y(:,:,2:nReads) = nan;
try
    for r = 2:nReads
        tiffObj.nextDirectory;
        Y(:,:,r) = tiffObj.read;
    end
catch ME %#ok<NASGU>
    readsuccess = false;
    return
end
Y = reshape(Y, size(Y,1), size(Y,2), nChannels, nFramesToRead);
readsuccess = true;
try
    tiffObj.nextDirectory;
catch
    disp('reached End of File');
    done = true;
end
end