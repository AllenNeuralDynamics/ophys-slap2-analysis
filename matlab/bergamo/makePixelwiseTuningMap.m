function makePixelwiseTuningMap(fns)

if nargin<1 || isempty(fns)
    fns = uipickfiles('FilterSpec','*_REGISTERED_RAW.tif');
end

if ~iscell(fns)
    fns = {fns};
end

for fix = 1:numel(fns)
    fn = fns{fix};
    disp(['Processing: ' fn])
    [dr,fn, ext] = fileparts(fn);
    fn = strcat(fn,ext);

    % get acquisition index
    [idxStart, idxEnd] = regexp(fn,'_\d*_');
    acqIdx = str2num(fn(idxStart+1:idxEnd-1));

    % get harp directory
    [idxStart, ~] = regexp(dr,'[\\/]scans[\\/]');
    harpDir = [dr(1:idxStart) 'harp'];

    photodiode = readmatrix(fullfile(harpDir,['photodiode_' num2str(acqIdx-1) '.csv']),'NumHeaderLines',1);
    orientations = readmatrix(fullfile(harpDir,['orientations_' num2str(acqIdx-1) '.csv']));
    
    opts = delimitedTextImportOptions("NumVariables", 2);
    
    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["Flyback", "Timestamp"];
    opts.VariableTypes = ["categorical", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, "Flyback", "EmptyFieldRule", "auto");
    
    % Import the data
    frameClk = readtable(fullfile(harpDir,['frameClk_' num2str(acqIdx-1) '.csv']), opts);
    
    clear numIdx opts
    
    frameClkSignal = table2array(frameClk(:,1))=='True';
    diffFrameClk = diff(frameClkSignal);
    
    frameClkTrue = frameClk([false; diffFrameClk == 1],2);
    frameClkTrue = table2array(frameClkTrue);
    disp([num2str(length(frameClkTrue)) ' flybacks detected']);
    
    uniqueOrientations = sort(unique(round(orientations * 360 / 2 / pi)));
    
    pd_mixture = fitgmdist(photodiode(:,2),2);
    [~, onCluster] = max(pd_mixture.mu);
    pdOn = cluster(pd_mixture,photodiode(:,2)) == onCluster;
    
    onTimes = photodiode(diff(pdOn) > 0,1);
    offTimes = photodiode(diff(pdOn) < 0,1);
    offTimes = [offTimes; photodiode(end,1)];

    % read alignment data
    [idxStart, ~] = regexp(fn,'_REGISTERED_RAW.tif');
    load(fullfile(dr,[fn(1:idxStart) 'ALIGNMENTDATA.mat']));

    [Ad, ~, ~] = networkScanImageTiffReader(fullfile(dr,fn));
    [nRow, nCol, ~] = size(Ad);
    disp([num2str(size(Ad,3)) ' frames read']);

    Ad = reshape(Ad, size(Ad,1), size(Ad,2), aData.numChannels, []);
    IM = reshape(Ad(:,:,end,:),[],size(Ad,4));
    clear Ad;

    th = graythresh((aData.recNegErr-min(aData.recNegErr)) ./ (max(aData.recNegErr)-min(aData.recNegErr))) .* (max(aData.recNegErr)-min(aData.recNegErr)) + min(aData.recNegErr);
    dsFac = floor(length(aData.motionC) / length(aData.motionDSc));
    recNegErrFull = interp1(1.5:dsFac:(1.5 + dsFac * (length(aData.recNegErr)-1)),aData.recNegErr,1:size(IM,2),'linear',median(aData.recNegErr));
    
    if mean(recNegErrFull > th) < 0.1
        IM(:,recNegErrFull > th) = nan;
    end

    validPix = find(mean(isnan(IM),2) < 0.2);
    IMsel = IM(validPix,:);

    IMsmooth = smoothdata(IMsel,2, 'movmean', 2*round(size(IMsel,2)/8)+1, 'omitnan');
    IMsel(isnan(IMsel)) = IMsmooth(isnan(IMsel));
    IMsel(isnan(IMsel)) = 0; %anything left over
    
    %%
    F0 = smoothExp(IMsel', 'movmedian', 10001)';
    
    %%
    dF = nan(size(IM));
    dF(validPix,:) = (IM(validPix,:) - F0);

    nanPix = mean(isnan(IM),2) > 0.5;

    pixResponse = nan([nRow*nCol, length(orientations) / length(uniqueOrientations), length(uniqueOrientations)]);

    for i = 1:length(orientations)
        startTime = onTimes(i);
        endTime = offTimes(i);
    
        trialIdxs = find((frameClkTrue > startTime) & (frameClkTrue < endTime));
        if i == 1
            offIdxs = find((frameClkTrue > startTime - 1) & (frameClkTrue < startTime));
        else
            offIdxs = find((frameClkTrue > offTimes(i-1)) & (frameClkTrue < startTime));
        end
        
        pixResponse(:, floor((i-1) / length(uniqueOrientations))+1, uniqueOrientations == round(orientations(i) * 360 / 2 / pi)) ...
            = max(0,mean(dF(:,trialIdxs),2,'omitnan') - mean(dF(:,offIdxs),2,'omitnan'));
        %mean(IM(:,trialIdxs),2,'omitnan') - mean(IM(:,offIdxs),2,'omitnan');

        pixResponse(mean(isnan(dF(:,trialIdxs)),2) > 0.1 | mean(isnan(dF(:,offIdxs)),2) > 0.1, ...
            floor((i-1) / length(uniqueOrientations))+1, uniqueOrientations == round(orientations(i) * 360 / 2 / pi)) ...
            = nan;
    end
    
    %%
    pixDir = nan(nRow*nCol,1);
    pixOSI = nan(nRow*nCol,1);
    
    pixAmp = sqrt(sum(squeeze(mean(pixResponse,2,'omitnan')).^2,2,'omitmissing'));
    
    for i = 1:size(pixResponse,1)
        tmp = [0; 0];
        for j = 1:size(pixResponse,3)
            tmp = tmp + squeeze(sum(pixResponse(i,:,j),'omitmissing')) .* [cosd(uniqueOrientations(j) * 2); sind(uniqueOrientations(j) * 2)];
        end
        
        pixOSI(i) = norm(tmp);
        if sum(abs(pixResponse(i,:,:)),'all','omitmissing') > 0; pixOSI(i) = pixOSI(i) ./ sum(abs(pixResponse(i,:,:)),'all','omitmissing'); end
        
        pixDir(i) = atan2d(tmp(2),tmp(1));
    end
    pixDir(pixDir < 0) = pixDir(pixDir < 0) + 360;
    
    
    pixDir(nanPix) = nan;
    pixOSI(nanPix) = nan;
    pixAmp(nanPix) = nan;
    
    
    imDir = reshape(pixDir,nRow,nCol);
    imOSI = reshape(pixOSI,nRow,nCol);
    imAmp = reshape(pixAmp,nRow,nCol);
    
    % Normalize the images to [0, 1]
    imDir = imDir / 360;
    
    osiMaxLim = prctile(pixOSI,99.9);
    imOSI = max(0, imOSI - min(imOSI(:)));
    imOSI = min(1, imOSI / (osiMaxLim - min(pixOSI)) );
    
    ampMaxLim = prctile(pixAmp,99.9);
    imAmp = max(0, imAmp - min(imAmp(:)));
    imAmp = min(1, imAmp / (ampMaxLim - min(pixAmp)) );
    
    % Create the HSV image
    hue = imDir; % Hue corresponds to the direction values
    saturation = imOSI; % Saturation corresponds to the OSI values
    value = imAmp; % Set value (brightness) to be equal to saturation
    
    % Combine into an HSV image
    hsvImage = cat(3, hue, saturation, value);
    
    % Convert HSV image to RGB
    rgbImage = hsv2rgb(hsvImage);
    
    % Create a figure with two subplots: one for the image and one for the color bar
    figure;
    
    % Display the image
    subplot(1, 2, 1);
    imshow(rgbImage);
    title('Pixelwise Tuning Map')
    
    % Generate the 2D color bar
    [hueGrid, satGrid] = meshgrid(linspace(0, 1, 256), linspace(0, 1, 256));
    valueGrid = ones(size(satGrid)); % Value (brightness) set equal to saturation
    
    % Combine into an HSV image for the color bar
    hsvColorBar = cat(3, hueGrid, satGrid, valueGrid);
    
    % Convert HSV color bar to RGB
    rgbColorBar = hsv2rgb(hsvColorBar);
    
    % Display the 2D color bar
    subplot(1, 2, 2);
    imshow(rgbColorBar);
    axis on;
    xlabel('Orientation');
    ylabel('OSI');
    title('2D Color Bar');
    
    % Adjust the color bar to match the hue representation
    xticks(linspace(0, 256, 5));
    xticklabels(arrayfun(@(x) sprintf('%d', x), linspace(0, 180, 5), 'UniformOutput', false));
    yticks(linspace(0, 256, 6));
    yticklabels(arrayfun(@(x) sprintf('%.1f', x), linspace(min(pixOSI), osiMaxLim, 6), 'UniformOutput', false));

    savefig(fullfile(dr,'pixelTuningMap.fig'));

end

end