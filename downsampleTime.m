function Y = downsampleTime(Y, ds_time)
    for ix = 1:ds_time
        Y = Y(:,:,:,1:2:(2*floor(end/2)))+ Y(:,:,:,2:2:end);
    end
end