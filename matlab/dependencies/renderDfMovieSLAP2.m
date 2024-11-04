function renderDfMovieSLAP2
%render nice movies from recordings of GluSnFR

dsfac = 2;
nChan = 2;
window_s = 0.25;


%load files
if ~nargin
    [fns,dr] = uigetfile('*DMD1*Hz*.tif', 'Select TIFF file', 'render Movie', 'MultiSelect','on') ;% '\\dm11\genie\Jihong_GENIE_stuff\iGluSnFR_spontaneous_activity' %'\\dm11\podgorskilab\GluSnFR4';
end
if ~iscell(fns)
    fns = {fns};
end

for fnIx = 1:length(fns)
    fn1 = fns{fnIx};
    IM1 = double(networkScanImageTiffReader([dr filesep fn1]));

    frexp = regexp(fn1, '[\d]+Hz', 'match');
    if isempty(frexp)
        error('Expected downsampled tiff filename to include the sampling rate as XXHz');
    end
    frameRate = str2double(frexp{1}(1:end-2));
    window = ceil(window_s*frameRate);

    fn2 = strrep(fn1,'DMD1', 'DMD2');
    IM2 = double(networkScanImageTiffReader([dr filesep fn2]));

    IM2(:, end+1:size(IM1,2),:) = nan;
    IM1(:, end+1:size(IM2,2),:) = nan;
    IM2(:, :,end+1:size(IM1,2)) = [];
    IM1(:, :, end+1:size(IM2,2)) = [];

    IM = cat(1,IM1,IM2); clear IM1 IM2
    IM = reshape(IM, size(IM,1), size(IM,2), nChan, []);

    sz = [size(IM,1) size(IM,2)];
    T = size(IM,4);
    [b1,a1] = butter(4, window/(T/2^dsfac)); %lowpass filter for defining F0

    for ch = [1]
        %downsample in time
        IMds = squeeze(IM(:,:,ch,:));
        for it = 1:dsfac
            IMds = IMds(:,:,1:2*floor(end/2));
            IMds=IMds(:,:,1:2:end) + IMds(:,:,2:2:end);
        end
        IMds = IMds./(2.^it); %keep scale;

        meanIM = mean(IMds,3, 'omitnan');
        
        % somaIM = ordfilt2(max(1,meanIM+1), [25 25]);
        % somaIM = imgaussfilt(somaIM, [15 15]);

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
        F0ds(:,:,1:ceil(window)) = repmat(F0ds(:,:,ceil(window)+1), 1, 1, ceil(window));
        clear e1;

        %compute dF
        dFds = IMds-F0ds;

        if ch==1
            %somaIM = repmat(sqrt(somaIM),1,1,T2);
            somaIM = max(50,meanIM+50);
            act = medfilt3(dFds, [3,3,1]);
            act = max(0, act - median(act(~nans)));
            act(repmat(setnan, 1,1,T2)) = 0;
            
            act = act./somaIM;
            dFmax = prctile(act(~nans), 99.97);
            act = act./dFmax;

            st = sqrt(max(0,F0ds./somaIM));
            st = st-prctile(st(~nans),10);
            st = 0.75*(st./prctile(st(~nans), 99.99));
            st(repmat(setnan, 1,1,T2)) = 0;

            st1 = st;
            act1 = act;
        else
            minSoma = median(meanIM(meanIM>0));
            somaIM = sqrt(max(0, meanIM-minSoma));
            
            %somaIM = repmat(sqrt(somaIM),1,1,T2);
            act = imgaussfilt(max(0, dFds), [3,3]);
            act(repmat(setnan, 1,1,T2)) = 0;
            
            act = act.*somaIM;
            dFmax = prctile(act(~nans), 99.9);
            act = act./dFmax;

            st = sqrt(max(0,F0ds./somaIM));
            st = st-prctile(st(:),1);
            st = 0.8*(st./prctile(st(:), 99.99));

            st2 = st;
            act2 = act;
        end
    end



    %plot flashes in red and structure in cyan (gamma corrected)
    %save to avi file
    frametime = (1/80)*(2.^dsfac);
    v = VideoWriter([dr filesep fn1(1:end-4) 'SUMMARY2.mp4'],'MPEG-4');
    v.Quality = 99;
    v.FrameRate = 10;
    open(v);
    hF = figure;
    hAx = axes;
    fix = 1;
    gray = st1(:,:,fix);
    red = act1(:,:,fix);
    F = max(0, min(1, cat(3, gray.*(1-red)+red, gray.*(1-red), gray.*(1-red))));
    hObj = imshow(F, 'InitialMagnification',500);
    hold on
    hT = text(10,10, ['t=' num2str((fix-1)*frametime, '%5.3f') ' s'], 'color', 'w');
    for fix = 1:size(st,3)
        gray = st1(:,:,fix);
        red = max(0,min(1,act1(:,:,fix)));
        cyan = max(0,min(1, st1(:,:,fix)));
        %F = max(0, min(1, cat(3, gray.*(1-red)+red, gray.*(1-red), gray.*(1-red))));
        F = max(0, min(1, cat(3, red, cyan.*sqrt(1-red), cyan.*sqrt(1-red))));
        hObj.CData = F;
        hT.String =  ['t=' num2str((fix-1)*frametime, '%5.3f') ' s'];
        fr = getframe(hAx);
        writeVideo(v,fr)
    end
    close(v)
    close(hF)

    clear st act mm st1 st2 act1 act2 IMds F0ds dFds v IM nans
end
end