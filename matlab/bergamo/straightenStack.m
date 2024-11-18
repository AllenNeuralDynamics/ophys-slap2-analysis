function [IMout]= straightenStack(IM, IMc, IMsk)
clipEdge = 40;
maxShift = 16;
scalingFac = round(size(IM,1)/size(IMc,1));
Y = makeHighPass(sum(IM(clipEdge+1:end-clipEdge,clipEdge+1:end-clipEdge,:,:),3)); % data order is XYCZ
Z2fft = fft2(Y(:,:,1,1));
IMout = IM;
dX = 0; dY = 0;
for Z1 = 1:(size(IM, 4)-1)
    Z1fft = Z2fft;
    Z2fft = fft2(Y(:,:,1,Z1+1));
    output = dftregistration_clipped(Z1fft,Z2fft,4,maxShift);
    dX = dX + output(4); dY = dY+output(3);
    IMout(:,:,:,Z1+1) = imtranslate(IM(:,:,:,Z1+1), [dX dY]);
    % if nargout>1
    %     IMc(:,:,:,Z1+1) = imtranslate(IMc(:,:,:,Z1+1), [dX dY]/scalingFac);
    % end
    % if nargout>2
    %     IMsk(:,:,:,Z1+1) = imtranslate(IMsk(:,:,:,Z1+1), [dX dY]/scalingFac);
    % end
end
end