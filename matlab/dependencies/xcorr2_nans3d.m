function [motion, R] = xcorr2_nans3d(frame, template3d, shiftsCenter, dShift)
%Kaspar Podgorski 2023
%perform a somewhat-efficient local normalized crosscorrelation for images with
%nans

%Adapted to use 3D template, Michael Xie 2024

%template: the template
%frame: the frame to be aligned; this has more NaNs

%shiftsCenter: the center offset around which to perform a local search
%dShift: the maximum shift (scalar, in pixels) to consider on each axis around shiftsCenter

dShift = round(dShift(1)); %sanity check

SE = strel(ones(2*dShift+1)); %we will erode with this structuring element to select the set of pixels that are always non-nan within the search space

fValid = ~isnan(frame) & circshift(~imdilate(any(isnan(template3d),3),SE), shiftsCenter); %valid pixels of the new frame
fValid(1:dShift,:) = false; fValid(end-dShift+1:end,:) = false; %remove edges
fValid(:,1:dShift) = false; fValid(:,end-dShift+1:end) = false; %remove edges

tValid = circshift(fValid, -shiftsCenter);

%fValid = imerode(~isnan(template) & circshift(~isnan(frame), -shiftsCenter), SE); %valid pixels of the fixed image

%mValid = circshift(fValid, shiftsCenter); %valid pixels of the moving image

F = frame(fValid); %fixed data;
ssF = sqrt(sum(F.^2));

%correlation is sum(A.*B)./(sqrt(ssA)*sqrt(ssB)); ssB is constant though
shifts =  -dShift:dShift;
C = nan(length(shifts), length(shifts), size(template3d,3));
% del = 5 * mad(F,1);
for drix = 1:length(shifts)
    for dcix = 1:length(shifts)
        for dzix = 1:size(template3d,3)
            T = template3d(:,:,dzix);
            T = T(circshift(tValid, -[shifts(drix) shifts(dcix)]));

            % a = abs(F-T);
            % clipVals = (a > del);
            % a(clipVals) = del * (a(clipVals) - del/2);
            % a(~clipVals) = a(~clipVals).^2/2;
            % 
            % C(drix,dcix) = -sum(a);
            
            ssT = sum(T.^2);
            C(drix,dcix,dzix) = sum(F .* T)./sqrt(ssT);
        end
    end
end

%find maximum of correlation map
[maxval,I] = max(C(:));
[rr,cc,zz] = ind2sub(size(C),I);

R = maxval./ssF; %correlation coefficient

if rr>1 && rr<length(shifts) && cc>1 && cc<length(shifts) && zz>1 && zz<size(template3d,3)
    %perform superresolution upsampling
    ratioR = min(1e6,(C(rr,cc,zz) - C(rr-1,cc,zz))/(C(rr,cc,zz) - C(rr+1,cc,zz)));
    dR = (1-ratioR)/(1+ratioR)/2;
    
    ratioC =min(1e6, (C(rr,cc,zz) - C(rr,cc-1,zz))/(C(rr,cc,zz) - C(rr,cc+1,zz)));
    dC = (1-ratioC)/(1+ratioC)/2;

    ratioZ =min(1e6, (C(rr,cc,zz) - C(rr,cc,zz-1))/(C(rr,cc,zz) - C(rr,cc,zz+1)));
    dZ = (1-ratioZ)/(1+ratioZ)/2;

    motion = [shiftsCenter' 0] + [shifts(rr)-dR shifts(cc)-dC zz-dZ];
else %the optimum is at an edge of search range; no superresolution
    motion = [shiftsCenter' 0] + [shifts(rr) shifts(cc) zz];
end

if any(isnan(motion))
    keyboard
end