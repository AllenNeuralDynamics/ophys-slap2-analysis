function renderDfMovie
%render nice movies from recordings of GluSnFR

dsfac = 2;
nChan = 2;
[b1,a1] = butter(4, 0.02); %lowpass filter for defining F0/bleaching correction


%load files
if ~nargin
    [fns,dr] = uigetfile('*.tif', 'Select TIFF file', 'render Movie', 'MultiSelect','on') ;% '\\dm11\genie\Jihong_GENIE_stuff\iGluSnFR_spontaneous_activity' %'\\dm11\podgorskilab\GluSnFR4';
end
if ~iscell(fns)
    fns = {fns}
end

for fnIx = 1:length(fns)
    fn = fns{fnIx};
    A = ScanImageTiffReader([dr filesep fn]);
    IM = double(A.data);

    IM = reshape(IM, size(IM,1), size(IM,2), nChan, []);
    clear f

    sz = [size(IM,1) size(IM,2)];
    T = size(IM,4);

    %downsample in time
    IMds = squeeze(IM(:,:,1,:));
    for it = 1:dsfac
        IMds = IMds(:,:,1:2*floor(end/2));
        IMds=IMds(:,:,1:2:end) + IMds(:,:,2:2:end);
    end
    IMds = IMds./(2.^it); %keep scale;

    meanIM = mean(IMds,3, 'omitnan');
    setnan = mean(isnan(IMds),3)>0.7;
    meanIM(setnan) = 0;
    T2 = size(IMds,3);

    %remove nans
    mm = repmat(meanIM, 1,1,T2);
    nans = isnan(IMds);
    IMds(nans) = mm(nans);

    %compute F0
    e1 = smoothdata(IMds,3, 'gaussian', 3);
    F0ds = medfilt3(e1, [1 1 25], 'replicate');
    for iter = 1:3 %make F0 smooth; this is slow
        F0ds = permute(filtfilt(b1,a1,permute(min(e1,F0ds), [3 1 2])), [2 3 1]);
    end
    clear e1;

    %compute F0 and dF
    dFds = IMds-F0ds;
    %F0 = reshape(interp1((0.5:1:T2).*(2.^it), reshape(F0ds,[],T2)',  1:T, 'linear', 'extrap')', size(IM));
    %dF = IM - F0;

    %ensure that F0 is always>0; We will additionally regularize when
    %calculating DFF
    % F0 = F0 - min(0, min(F0,[],3));
    % clear F0ds IM
    st = max(0,F0ds);
    act = max(0, dFds - median(dFds(~nans)));
    act(repmat(setnan, 1,1,T2)) = 0;
    st = sqrt(st./prctile(st(:), 99.99));
    st = st-prctile(st(:),1);

    dFmax = prctile(act(:), 99.95);
    act = act./dFmax;

    %plot flashes in red and structure in cyan (gamma corrected)
    %save to avi file
    v = VideoWriter([dr filesep fn(1:end-4) '.avi'],'Motion JPEG AVI');
    v.Quality = 95;
    v.FrameRate = 15;
    open(v);

    for fix = 1:size(st,3)
        F = max(0, min(1, cat(3, act(:,:,fix), st(:,:,fix), st(:,:,fix))));
        writeVideo(v,F)
    end
    close(v)
end