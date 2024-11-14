function assessResiduals(IMrec, IMsel, selPix2D)
nFrames = size(IMsel,2);

%for each source, compute the movie
mRec = nan([numel(selPix2D) size(IMsel,2)]);
mRec(selPix2D(:),:) = IMrec;
mRec = reshape(mRec, [size(selPix2D) nFrames]);

%reshape IMsel into the movie
mObs = nan([numel(selPix2D) size(IMsel,2)]);
mObs(selPix2D(:),:) = IMsel;
mObs = reshape(mObs, [size(selPix2D) nFrames]);

mObs = downsampleTime(mObs,2);
mRec = downsampleTime(mRec,2);

figure('name', 'Observed'), imshow3D(mObs);
figure('name', 'Reconstruction'), imshow3D(mRec)
figure('name', 'Residual'), imshow3D(mObs-mRec)

end

function Y = downsampleTime(Y, ds_time)
for ix = 1:ds_time
    Y = Y(:,:,1:2:(2*floor(end/2)))+ Y(:,:,2:2:end);
end
end
