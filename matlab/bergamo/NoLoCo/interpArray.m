function IMsel = interpArray (IM, sel, shiftRC)
%linearly interpolate the 3D matrix IM in each 2D plane at the selected
%pixels sel, shifted by shiftRC
%returns a 2D array: [sum(sel) x size(IM,3)]
sz = size(IM);
IMsel = nan(sum(sel(:)), sz(3));
IM = reshape(IM, sz(1)*sz(2), []);

inds = zeros(size(sel));
inds(sel) = 1:sum(sel(:));

sel = sel(1:sz(1), 1:sz(2));
inds = inds(1:sz(1), 1:sz(2));

intShift = floor(shiftRC);
shiftRC = shiftRC-intShift;
sel = imtranslate(sel, [intShift(2) intShift(1)]);
inds = imtranslate(inds, [intShift(2) intShift(1)]);

%ensure that all the masks have the same number of values:
if shiftRC(1)>0.05
    sel(end,:) = false; inds(end,:) = false;
end
if shiftRC(2)>0.05
    sel(:,end) = false; inds(:,end) = false;
end
inds = inds(sel);

mask00 = sel;                    %unshifted
mask10 = imtranslate(sel,[0 1]); %shifted 1 row
mask01 = imtranslate(sel,[1 0]); %shifted 1 col
mask11 = imtranslate(sel,[1 1]); %shifted 1 row and 1 col

if shiftRC(1)>0.05 %the subpixel shift is nonnegligible, so interpolate
    R0 = (1-shiftRC(1)).*IM(mask00(:),:) + shiftRC(1).*IM(mask10(:),:);
    R1 = (1-shiftRC(1)).*IM(mask01(:),:) + shiftRC(1).*IM(mask11(:),:);
else %subpixel shift is negligible, use the unshifted data (this prevents NaNing out good data at edges)
    R0 = IM(mask00(:),:);
    R1 = IM(mask01(:),:);
end
if shiftRC(2)>0.05
    IMsel(inds,:) = (1-shiftRC(2)).*R0 + shiftRC(2).*R1;
else
    IMsel(inds,:) = R0;
end
end