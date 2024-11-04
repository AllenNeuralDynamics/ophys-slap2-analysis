function  [rawIMs, meanIM, IMc, aData, peaks, discardFrames]= loadAndProcessTrialAsync(dr, fn, numChannels, params)
    
    nInitFrames = ceil(params.discardInitial_s/params.frametime); %initial frames to discard

    %load the tiff
    IM = copyReadDeleteScanImageTiff([dr filesep fn]);
    IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []); %deinterleave;

    meanIM= mean(IM,4, 'omitnan');
    nanPx = repmat(mean(isnan(IM), [3 4])>0.15, 1,1,size(IM,3));
    meanIM(nanPx) = nan;

    rawIMs = squeeze(IM(:,:,1,:));
    clear IM

    %load alignment data
    fnStemEnd = strfind(fn, '_REGISTERED') -1;
    load([dr filesep fn(1:fnStemEnd) '_ALIGNMENTDATA.mat'], 'aData');
    aData.dsFac = 1; %SLAP2 data does not have downsampling per se
    params.dsFac = aData.dsFac;
    params.frametime = aData.frametime;

    %discard motion frames
    tmp1 = aData.aRankCorrDS(:)-smoothExp(aData.aRankCorrDS(:),'movmedian', ceil(2/(aData.frametime*aData.dsFac))); %-smoothdata(aData.aRankCorrDS,2, 'movmedian', ceil(2/aData.frametime));
    tmp2 = aData.recNegErr(:)-smoothExp(aData.recNegErr(:),'movmedian', ceil(2/(aData.frametime*aData.dsFac))); %-smoothdata(aData.aRankCorrDS,2, 'movmedian', ceil(2/aData.frametime));
    tmp3 = aData.aError(:) -smoothExp(aData.aError(:),'movmedian', ceil(2/(aData.frametime*aData.dsFac)));

    tmp1 = tmp1./std(tmp1(nInitFrames+1:end)); 
    tmp2 = -tmp2./std(tmp2(nInitFrames+1:end)); 
    tmp3 = -tmp3./std(tmp3(nInitFrames+1:end)); 
    
    tmp = (tmp1 + tmp2 + tmp3)/6; %decrease denominator to be more stringent on motion correction
    discardFrames = imclose(imdilate(tmp<-2, ones(5,1)) | (tmp<-1 & imdilate(tmp<-2, ones(21,1))), ones(11,1));
    discardFrames(1:nInitFrames) = true;
    rawIMs(:,:,discardFrames) = nan;
    
    
    [IMc, peaks] = localizeSourcesSLAP2(rawIMs, aData, params);
end