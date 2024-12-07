function motionLUT = makeLookupTable(reference, sparseMaskInds, numFastZs, fastZ2RefZ, varargin)

[dmdPixelsPerColumn, dmdPixelsPerRow, numPlanes, numChannels] = size(reference);

if nargin == 3
    yMotRange = -25:25;
    xMotRange = -25:25;
    zMotRange = -10:10;
    channels = 1:numChannels;
else
    yMotRange = varargin{1};
    xMotRange = varargin{2};
    zMotRange = varargin{3};
    channels = varargin{4};
end

numY = length(yMotRange);
numX = length(xMotRange);
numZ = length(zMotRange);


fprintf("Calculating lookup table... ")

tic;
motionLUT = single(zeros(length(yMotRange), length(xMotRange), length(zMotRange),max(channels),max(sparseMaskInds(:,2)))); % X x Y x Z x channels x no. of superpixels
motionLUT(:) = nan;

refPixs = zeros(max(sparseMaskInds(:,2)),1);

parfor n = 1:max(sparseMaskInds(:,2))
    openPixs = sparseMaskInds(sparseMaskInds(:,2) == n, 1);
    refPix = ceil(length(openPixs)/2);
    refPixs(n) = openPixs(refPix);
end

[r, c, d] = ind2sub([dmdPixelsPerColumn dmdPixelsPerRow numFastZs], refPixs);
d = fastZ2RefZ(d);

for y = 1:numY
    for x = 1:numX
        for z = 1:numZ
            shiftedR = r + yMotRange(y);
            shiftedC = c + xMotRange(x);
            shiftedD = d + zMotRange(z);

            validInds = (shiftedR >= 1) & (shiftedR <= dmdPixelsPerColumn) & (shiftedC >= 1) & (shiftedC <= dmdPixelsPerRow) ...
                & (shiftedD >= 1) & (shiftedD <= numPlanes);

            if sum(validInds) == 0; continue; end

            % select expected value of superpixel at (x,y,z)
            % displacement from the reference image
            for chIx = channels
                motionLUT(y,x,z,chIx,validInds) = reference(sub2ind([dmdPixelsPerColumn dmdPixelsPerRow numPlanes max(channels)], ...
                    shiftedR(validInds), shiftedC(validInds), shiftedD(validInds),chIx * ones(size(validInds))));
            end
        end
    end
end

fprintf('done - took %f sec\n', toc);