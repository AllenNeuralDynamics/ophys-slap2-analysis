function [correlationTable, scalingFactorTable] = correlationLikelihoodTable(data, likelihood_means, ySearch, xSearch, zSearch, channels)

validSPs = ~any(isnan(data),2);
sub_likelihood_means = likelihood_means(ySearch,xSearch,zSearch,channels,:);

scalingFactorTable = ones(size(sub_likelihood_means,1:4));

for chIx = 1:length(channels)
    expectedData = squeeze(sub_likelihood_means(:,:,:,chIx,validSPs)); % [y, x, z, validSPs]
    tmpData = reshape(data(validSPs, chIx), [1, 1, 1, sum(validSPs)]);

    numerator = sum(tmpData .* ~isnan(expectedData), 4);
    denominator = sum(expectedData, 4, 'omitnan');
    scalingFactorTable(:,:,:,chIx) = numerator ./ denominator;
end

inds = find(validSPs);
nCh = length(channels);
slice = sub_likelihood_means(:,:,:, :, inds); % [ny nx nz nCh nValid]
[nY, nX, nZ, ~, nValid] = size(slice);

% Same ordering as exp_mat(:) for exp_mat = nCh x nValid: column-major on [ch x valid]
obs_flat = reshape(data(inds, :)', [], 1);
obs_c = obs_flat - mean(obs_flat, 'omitnan');
norm_obs = sqrt(sum(obs_c.^2, 'omitnan'));
if norm_obs == 0 || ~isfinite(norm_obs)
    correlationTable = nan(nY, nX, nZ);
else
    S = reshape(permute(slice, [4 5 1 2 3]), nCh * nValid, []);
    M = S.'; % nY*nX*nZ x (nCh*nValid), one joint expected vector per shift
    rowMean = mean(M, 2, 'omitnan');
    M_c = M - rowMean;
    dot_vec = M_c * obs_c;
    norm_exp = sqrt(sum(M_c.^2, 2, 'omitnan'));
    correlationTable = reshape(dot_vec ./ (norm_obs .* norm_exp + eps), nY, nX, nZ);
end

correlationTable(mean(isnan(sub_likelihood_means),[4 5]) == 1) = nan;
correlationTable(~isfinite(correlationTable)) = -1e10;

scalingFactorTable(~isfinite(scalingFactorTable)) = 1;

end
