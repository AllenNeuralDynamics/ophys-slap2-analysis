classdef spineAnalysis < handle
    properties
        hFigPanel; %image panel
        hButtonPanel; %button panel
        hF; %figure handle
        hGrid;
        hAx; %axes handle
        hButton1; hButton1c;
        hButton2;
        hButton3;
        hButton4;
        hButtonCh;
        hEditCh;

        hROIs;
        hIm; %image object

        chDim = 3; %images are in XYCT order
        timeDim = 4; %images are in XYCT order
        showChannel = 1;

        %butterworth filter parameters for highpass
        b2;
        a2;

        IM;
        IMc; %correlation image; normalized to median 0 standard deviation 1
        IMsk; %skewness image; normalized to median 0 standard deviation 1
        IMact; %the image used to display activity
        IMavg; %average image, normalized to [0 1]
        bleach;
        aData; %alignment data
        dsFac; %downsampling factor of the image file, used to properly set the framerate;
        numChannels;
        fn;
        dr = [];
        drsave;
        fnsave;
        fnStem;

        defaultChLabels = {'iGluSnFr'; 'jRGECO'};
        Clevel = 3; %adjustable threshold for display of correlation image, in standard deviations
    end

    methods
        function obj = spineAnalysis(arg1)
            if ~nargin || isempty(arg1)
                [obj.fn, obj.dr] = uigetfile('*DOWNSAMPLED*.tif');
            elseif isstring(arg1) || ischar(arg1)
                [obj.dr, fn] = fileparts(arg1);
                obj.fn = [fn '.tif'];
            else
                error('incorrect input!')
            end

            %load motion correction data
            fnStemEnd = strfind(obj.fn, '_REGISTERED') -1;
            obj.fnStem = obj.fn(1:fnStemEnd);
            try
                load([obj.dr filesep obj.fn(1:fnStemEnd) '_ALIGNMENTDATA.mat'], 'aData');
                obj.aData = aData;
                if isfield(aData, 'numChannels')
                    obj.numChannels = aData.numChannels;
                else
                    obj.numChannels = 1;
                end
            catch
                msgbox('Could not load alignment data; traces will not be decorrelated to motion.')
                obj.numChannels = 1;
            end

            %load image data
            A = ScanImageTiffReader([obj.dr filesep obj.fn]);
            IM = A.data;
            IM = reshape(IM, size(IM,1), size(IM,2), obj.numChannels, []); %deinterleave;
            obj.IM = permute(IM, [2 1 3 4]); % X Y C T order

            clear IM;

            %infer downsampling from filename
            str = extractBetween(obj.fn, '_DOWNSAMPLED-', 'x.tif');
            if ~isempty(str)
                obj.dsFac = str2double(str{1});
            else
                obj.dsFac = 1;
            end

            %initialize butterworth filter
            [obj.b2,obj.a2] = butter(4, 0.05, 'high'); %highpass filter for generating correlation image

            %compute the bleaching curve
            selValid = ~any(isnan(obj.IM(:,:,1,:)),obj.timeDim);
            tmp = reshape(obj.IM, [], obj.numChannels, size(obj.IM,obj.timeDim));
            offset = prctile(mean(tmp,3, 'omitnan'), 1, 1);
            obj.IM = obj.IM-reshape(offset, [1 1 numel(offset) 1]);
            tmp = tmp - offset;
            obj.bleach = mean(tmp(selValid,:,:,:), 1, 'omitnan');

            %compute average image
            IMavg = mean(obj.IM,obj.timeDim, 'omitnan');
            IMgamma = sqrt(max(0, IMavg./prctile(IMavg(:), 99.99)));
            obj.IMavg = IMgamma;

            %compute correlation image
            disp('computing correlation image');

            %select valid frames for cross-correlation analysis
            if ~isempty(obj.aData)
                err = aData.aError(1:obj.dsFac:end);
                errStd = sqrt(estimatenoise(err));
                valid = err < (ordfilt2(err, 40, ones(1,200), 'symmetric')+5*errStd);
            else
                valid = true(1, size(obj.IM,3));
            end


            noNan = double(obj.IM(:,:,:,valid));
            noNan = noNan - mean(noNan,obj.timeDim, 'omitnan');
            nanInds = isnan(noNan);
            noNan(nanInds)=0;
            HP = permute(filtfilt(obj.b2,obj.a2,permute(noNan, [obj.timeDim 1 2 obj.chDim])), [2 3 4 1]); clear noNan;
            HP(nanInds) = nan;

            %discard pixels with few measurements
            nanOut = sum(~isnan(HP(:,:,1,:)),4)<(size(HP,4)/3);
            HP(repmat(nanOut, [1 1 size(HP,3) size(HP,4)])) = nan;

            ss = sum(HP.^2,obj.timeDim, 'omitnan');
            vertC = sum(HP .* circshift(HP, [1 0 0 0]),obj.timeDim, 'omitnan')./sqrt(ss.*circshift(ss, [1 0 0 0]));
            %             vertN = sum(~(nanInds | circshift(nanInds, [1 0 0 0])), obj.timeDim);
            %             vertC = (vertC-median(vertC, [1 2], 'omitnan')).*sqrt(vertN);

            horzC = sum(HP .* circshift(HP, [0 1 0 0]),obj.timeDim, 'omitnan')./sqrt(ss.*circshift(ss, [0 1 0 0]));
            %             horzN = sum(~(nanInds | circshift(nanInds, [0 1 0 0])), obj.timeDim);
            %             horzC = (horzC - median(horzC, [1 2], 'omitnan')).*sqrt(horzN);

            C = mean(cat(obj.timeDim, horzC, circshift(horzC,1,2), vertC, circshift(vertC, 1,1)),obj.timeDim, 'omitnan');
            obj.IMc = (C-median(C, [1 2], 'omitnan'))./std(C, 0,[1 2],'omitnan');
            
            HP = imgaussfilt(HP, [0.6 0.6]);
            sk = skewness(HP,0,obj.timeDim).*IMgamma;
            obj.IMsk = (sk-median(sk(:), 'omitnan'))./std(sk, 0,'all','omitnan');
            obj.IMact = obj.IMc + obj.IMsk;

            %gui creation
            obj.hF = uifigure('name', ['Spine analysis - ' obj.fn], 'KeyPressFcn', @obj.kpf);
            obj.hGrid = uigridlayout(obj.hF);
            obj.hGrid.RowHeight = {50, '1x'};
            obj.hGrid.ColumnWidth = {'1x'};

            obj.hButtonPanel = uigridlayout(obj.hGrid);
            obj.hButtonPanel.RowHeight = {'1x'};
            obj.hButtonPanel.ColumnWidth = {100,100,100,100,100,100, 100,100,100};
            obj.hButton1 = uibutton(obj.hButtonPanel, 'Text', 'Draw Polygon', 'ButtonPushedFcn', @obj.drawPolygon);
            obj.hButton1c = uibutton(obj.hButtonPanel, 'Text', 'Draw Circle', 'ButtonPushedFcn', @obj.drawCircle);
            obj.hButton2 = uibutton(obj.hButtonPanel, 'Text', 'Save ROIs', 'ButtonPushedFcn', @obj.saveROIs);
            obj.hButton3 = uibutton(obj.hButtonPanel, 'Text', 'Load ROIs', 'ButtonPushedFcn', @obj.loadROIs);
            obj.hButton4 = uibutton(obj.hButtonPanel, 'Text', 'Analyze', 'ButtonPushedFcn', @obj.analyze);
            for chIx = 1:obj.numChannels
                obj.hButtonCh(chIx) = uibutton(obj.hButtonPanel, 'Text', ['Show Ch' num2str(chIx) ':'], 'ButtonPushedFcn', @(x1,x2)(obj.setCh(chIx)));
                obj.hEditCh(chIx) = uieditfield(obj.hButtonPanel, 'Value', obj.defaultChLabels{chIx});
            end

            obj.hAx = uiaxes(obj.hGrid);
            obj.hIm = imshow(cat(3, obj.IMact(:,:,obj.showChannel)./obj.Clevel, obj.IMavg(:,:,obj.showChannel), obj.IMavg(:,:,obj.showChannel)),'parent', obj.hAx);
            set(obj.hIm, 'HitTest', 'off', 'pickableparts', 'none', 'clipping', 'on')
            set(obj.hAx, 'hittest', 'off', 'PickableParts', 'all')
        end

        function setCh(obj,chIx)
            obj.showChannel = chIx;
            updateDisplay(obj);
        end

        function updateDisplay(obj)
            set(obj.hIm, 'CData',cat(3, obj.IMact(:,:,obj.showChannel)./obj.Clevel, obj.IMavg(:,:,obj.showChannel), obj.IMavg(:,:,obj.showChannel)))
        end

        function sData = analyze(obj, arg1, arg2)
            if ~isempty(obj.aData) %populate motion predictors
                motPreds = [obj.aData.motionC ; obj.aData.motionR ; obj.aData.motionC.^2 ; obj.aData.motionR.^2 ; obj.aData.motionC.^3 ; obj.aData.motionR.^3 ; obj.aData.motionC.*obj.aData.motionR ; ones(size(obj.aData.motionC))];
            end

            %generate masks for the ROIs
            tix = 0;
            sData.names ={}; %INITIALIZE
            for rix = 1:length(obj.hROIs)
                try
                    roiObj = findobj(obj.hROIs(rix));
                catch
                    continue
                end

                tix = tix+1; %tix will count the total # of valid ROIs
                name = roiObj.Label;
                assert(~isempty(name), 'ROI names cannot be empty!');
                assert(~any(strcmpi(sData.names, name)), ['A ROI name was duplicated! ROI: ' name])
                sData.names{tix} = name;

                sData.mask{tix} = roiObj.createMask;
                [rr,cc] = find(sData.mask{tix});
                sData.bbox{tix} = [min(cc) min(rr) ; max(cc) max(rr)];
                Dmask{tix} = sData.mask{tix}(sData.bbox{tix}(1,2):sData.bbox{tix}(2,2), sData.bbox{tix}(1,1):sData.bbox{tix}(2,1)); %#ok<AGROW> 
            end
            nROIs = tix;

            %load the high time resolution data
            ind =strfind(obj.fn, '_DOWNSAMPLED');
            if ~isempty(ind)
                disp('Reading High-res data')
                fnRaw = [obj.fn(1:ind) 'RAW.tif'];
                nframes = obj.dsFac * size(obj.IM,obj.timeDim); %/ obj.numChannels;

                %initialize data series
                for rix =nROIs:-1:1
                    D{rix} = nan([sData.bbox{rix}(2,2)-sData.bbox{rix}(1,2)+1, sData.bbox{rix}(2,1)-sData.bbox{rix}(1,1)+1, obj.numChannels, nframes]);
                end

                finfo = dir([obj.dr filesep fnRaw]);
                if finfo(1).bytes<1e10
                    %file is small enough to fit into memory
                    A = ScanImageTiffReader([obj.dr filesep fnRaw]);
                    fdata = permute(A.data(), [2 1 3]);
                    fdata = reshape(fdata, size(fdata,1),size(fdata,2), obj.numChannels, nframes);
                    for rix = 1:nROIs
                        D{rix} = fdata(sData.bbox{rix}(1,2):sData.bbox{rix}(2,2), sData.bbox{rix}(1,1):sData.bbox{rix}(2,1),:,:);
                    end
                elseif nframes<2^16
                    t = Tiff([obj.dr filesep fnRaw], 'r');
                    t.setDirectory(1);
                    try
                        for fix = 1:nframes
                            for cix = 1:obj.numChannels
                                %pageInd = (fix-1)*obj.numChannels + cix
                                %fdata = imread([obj.dr filesep fnRaw], pageInd);
                                fdata = t.read;
                                for rix = 1:nROIs
                                    D{rix}(:,:,cix,fix) = fdata(sData.bbox{rix}(1,2):sData.bbox{rix}(2,2), sData.bbox{rix}(1,1):sData.bbox{rix}(2,1));
                                end
                                t.nextDirectory;
                            end
                        end
                    catch ME %#ok<NASGU>
                        disp(['Read ' int2str(fix) ' frames'])
                    end
                else
                    %file is too big to fit into memory and has too many
                    %pages for Tiff library
                    disp('file is too big to fit into memory, and too many pages for Tiff library... will solve this later... trying anyway')
                    t = Tiff([obj.dr filesep fnRaw], 'r');
                    t.setDirectory(1);
                    try
                        for fix = 1:nframes
                            for cix = 1:obj.numChannels
                                %pageInd = (fix-1)*obj.numChannels + cix
                                %fdata = imread([obj.dr filesep fnRaw], pageInd);
                                fdata = t.read;
                                for rix = 1:nROIs
                                    D{rix}(:,:,cix,fix) = fdata(sData.bbox{rix}(1,2):sData.bbox{rix}(2,2), sData.bbox{rix}(1,1):sData.bbox{rix}(2,1));
                                end
                                t.nextDirectory;
                            end
                        end
                    catch ME %#ok<NASGU>
                        disp(['Read ' int2str(fix) ' frames'])
                    end
                end
            end

            sData.traceAvg = nan(nframes, nROIs, obj.numChannels);
            sData.trace0 = nan(nframes, nROIs, obj.numChannels);
            sData.trace1 = nan(nframes, nROIs, obj.numChannels);
            sData.trace2 = nan(nframes, nROIs, obj.numChannels);
            for rix = 1:nROIs
                Dsz = size(D{rix});
                for cix = 1:obj.numChannels

                    DD = double(reshape(D{rix}(:,:,cix,:), [], Dsz(4)));
                    nans = isnan(DD);
                    numVal = sum(~nans, 2);
                    [~,p] = sort(numVal,'descend');
                    r = (1:length(p))'; r(p) = r;

                    %compute average time trace within the ROI, ignoring nans
                    traceAvg = mean(DD,1, 'omitnan');

                    %select pixels and frames to include
                    selpix = (numVal./max(numVal))>(r./max(r) + 0.05);
                    selframes = ~any(nans(selpix,:),1);
                    selpix = ~any(nans(:,selframes),2);
                    DD = DD(selpix,selframes); %discard nans

                    trace0= sum(DD(Dmask{rix}(selpix),:),1); %trace0 is the SUM over all selected pixels

                    %remove motion noise in the singular components of the movie
                    mDD = mean(DD,2);
                    [U,S,V] = svds(DD-mDD, floor(sum(selpix)/10));
                    Vcorr = V;
                    for vix = 1:size(V,2)
                        [~, Vcorr(:, vix)] = obj.decorrelateMotion(V(:,vix), motPreds(:,selframes));
                    end
                    DDcorr = DD-(U*S*Vcorr'); %motion-corrected data: the original data minus the reconstructed motion artefact
                    trace1= mean(DDcorr(Dmask{rix}(selpix), :),1); % the trace computed from the mask alone

                    %perform NMF
                    nfacs = max(2, floor(sum(selpix)/12));
                    try
                        [W,H] = nnmf(double(DDcorr),nfacs);

                        %find NMF components spatially correlated to trace0, temporally correlated, and with positive skewness
                        HP = filtfilt(obj.b2,obj.a2,H')';
                        t0HP = filtfilt(obj.b2,obj.a2,trace0);
                        Cspace = corr(W, Dmask{rix}(selpix));
                        Ctime = corr(HP', t0HP'); %high frequency correlation
                        SK = skewness(HP, 0, 2);
                        sel = Cspace>0.1 & Ctime>0.1 & SK>0.2;
                        trace2 = mean(W(:,sel)*H(sel,:),1);
                    catch
                        trace2 = NaN(size(trace1));
                    end

                    %compile output
                    sData.traceAvg(:,rix,cix) = traceAvg; %average value within ROI
                    sData.trace0(selframes, rix,cix) = trace0; %raw ROI, without decorrelating motion variables
                    sData.n0(rix,cix) = estimatenoise(trace0); %noise for trace0
                    sData.trace1(selframes, rix,cix) = trace1; %raw ROI, with motion variables decorrelated
                    sData.n1(rix,cix) = estimatenoise(trace1); %noise for trace1
                    sData.trace2(selframes, rix,cix) = trace2; %NMF-based traces
                    sData.n2(rix,cix) = estimatenoise(trace2); %noise for trace2
                    sData.contamination(rix, cix) = 1-corr(trace0', trace2'); %a measure of how contaminated the raw trace is by other factors/noise
                end
            end

            sData.frametime = obj.aData.frametime;
            sData.dsFac = obj.dsFac;
            sData.bleach = obj.bleach; %bleaching curve for this field of view, in downsampled time

            %plot the data
            colors = [0.8 0.3 0.3; 0.25 0.5 0.25; 0.25 0.25 0.5];
            T = obj.aData.frametime * (1:size(sData.trace0,1)); %time vector
            figure,
            for pix = 1:nROIs
                for cix = 1:obj.numChannels
                    hAx(pix,cix) = subplot(nROIs*obj.numChannels, 1, (pix-1)*obj.numChannels+cix);
                    tNorm = sData.trace1(:,pix,cix)./sqrt(sData.n1(rix,cix)); %normalized trace
                    plot(T, tNorm, 'color', sqrt(colors(cix,:))); hold on;
                    plot(T, nanfastsmooth(tNorm, 7, 3, 0.5), 'color', colors(cix,:), 'linewidth', 2);
                    ylabel([sData.names{pix} ' - Ch' int2str(cix)]);
                    if pix<nROIs
                        set(hAx(pix,cix), 'xticklabel', [], 'tickdir', 'out');
                    end
                end
            end
            xlabel('time (s)');
            linkaxes(hAx, 'x');

            %save output
            if isempty(obj.drsave) || ~any(obj.drsave)
                obj.drsave = obj.dr;
            end
            % [obj.fnsave, obj.drsave] = uiputfile([obj.drsave filesep obj.fnStem '_TRACES.h5']);
            obj.fnsave = [obj.fnStem '_TRACES.h5'];
            
            if obj.fnsave %if user didn't cancel
                obj.saveAsH5([obj.drsave filesep obj.fnsave(1:end-3) '.h5'] , sData);
            end

            %save figure
            disp('saving figure')
            exportgraphics(obj.hAx, [obj.drsave filesep obj.fnStem '_FIGURE.pdf'], 'ContentType', 'vector');

        end

        function saveAsH5(obj, fname, sData)
            if exist(fname, 'file')
                resp = questdlg(['The file ' fname ' already exists. Overwrite?'], 'Overwrite?', 'Yes', 'No', 'No');
                switch resp
                    case 'Yes'
                        delete(fname);
                    case 'No'
                        return
                end
            end
            nROIs = length(sData.names);

            szTrace0 = [size(sData.trace0,1) obj.numChannels]; %size of each trace
            szTrace1 = [size(sData.trace1,1) obj.numChannels];


            %saves spine analyses data in H5 format for easier python analysis
            for roiIx = 1:nROIs
                roiName = sData.names{roiIx};

                %create data
                h5create(fname,['/fluo/raw/' roiName], [szTrace0(2) szTrace0(1)], 'Datatype', 'single', 'ChunkSize', min([szTrace0(2) szTrace0(1)], [1 2000]), 'Deflate', 5); % fluorescence from this ROI
                h5create(fname,['/fluo/mcorr/' roiName], [szTrace1(2) szTrace1(1)], 'Datatype', 'single', 'ChunkSize', min([szTrace1(2) szTrace1(1)], [1 2000]), 'Deflate', 5); % fluorescence from this ROI
                h5create(fname,['/fluo/avgwNaN/' roiName], [szTrace1(2) szTrace1(1)], 'Datatype', 'single', 'ChunkSize', min([szTrace1(2) szTrace1(1)], [1 2000]), 'Deflate', 5); % fluorescence from this ROI

                % write data
                h5write(fname,['/fluo/raw/' roiName],permute(sData.trace0(:,roiIx,:), [3 1 2]));
                h5write(fname,['/fluo/mcorr/' roiName],permute(sData.trace1(:,roiIx,:), [3 1 2]));
                h5write(fname,['/fluo/avgwNaN/' roiName],permute(sData.traceAvg(:,roiIx,:), [3 1 2]));

                %write attributes
                h5writeatt(fname,['/fluo/raw/' roiName],'fs', 1/obj.aData.frametime); %sampling frequency
                h5writeatt(fname,['/fluo/raw/' roiName],'noise', sData.n0(roiIx,:)); %estimated noise, per channel
                h5writeatt(fname,['/fluo/mcorr/' roiName],'fs', 1/obj.aData.frametime); %sampling frequency
                h5writeatt(fname,['/fluo/mcorr/' roiName],'noise', sData.n1(roiIx,:)); %estimated noise, per channel
                h5writeatt(fname,['/fluo/avgwNaN/' roiName],'fs', 1/obj.aData.frametime); %sampling frequency

                chans = {};
                for chIx = 1:obj.numChannels
                    chans = cat(1, chans, {get(obj.hEditCh(chIx), 'Value')});
                end
                h5writeatt(fname,['/fluo/raw/' roiName],'chans', chans); %channel names
                h5writeatt(fname,['/fluo/raw/' roiName],'nPix', sum(sData.mask{roiIx}(:))); %channel names
                h5writeatt(fname,['/fluo/mcorr/' roiName],'chans', chans); %channel names
                h5writeatt(fname,['/fluo/mcorr/' roiName],'nPix', sum(sData.mask{roiIx}(:))); %channel names
                h5writeatt(fname,['/fluo/avgwNaN/' roiName],'chans', chans); %channel names
            end

            %save global information
            h5create(fname,"/fluo/motionC",size(obj.aData.motionC), 'Datatype', 'single','ChunkSize', size(obj.aData.motionC), 'Deflate', 5); % XY movement
            h5write(fname,"/fluo/motionC",obj.aData.motionC);
            h5create(fname,"/fluo/motionR",size(obj.aData.motionR), 'Datatype', 'single','ChunkSize', size(obj.aData.motionR), 'Deflate', 5); % XY movement
            h5write(fname,"/fluo/motionC",obj.aData.motionR);
            h5create(fname,"/fluo/alignmentError",size(obj.aData.aError), 'Datatype', 'single','ChunkSize', size(obj.aData.aError), 'Deflate', 5); % XY movement
            h5write(fname,"/fluo/alignmentError",obj.aData.aError);
            
            h5writeatt(fname,"/fluo/motionR",'fs', 1/obj.aData.frametime);
            h5writeatt(fname,"/fluo/motionC",'fs', 1/obj.aData.frametime);
            h5writeatt(fname,"/fluo/alignmentError",'fs', 1/obj.aData.frametime);

        end

        function saveROIs (obj, arg1, arg2)
            if isempty(obj.drsave)
                obj.drsave = obj.dr;
            end
            [obj.fnsave, obj.drsave] = uiputfile([obj.drsave filesep obj.fnStem '_ROIs.mat']);
            if isempty(obj.fnsave) || ~any(obj.fnsave)
                return
            end
            roiData = {};
            for rix = 1:length(obj.hROIs)
                try
                    S = rmfield(get(obj.hROIs(rix)), {'Parent', 'Children', 'ContextMenu'});
                    roiData = [roiData {S}];
                catch
                    %this object didn't exist
                end
            end
            save([obj.drsave filesep obj.fnsave], 'roiData');
        end

        function [corrected, traceMot] = decorrelateMotion(obj, trace, motPreds)

            %highpass filter trace
            trace = trace(:) - nanfastsmooth(trace(:), min(length(trace)/2, ceil(2/obj.aData.frametime)), 3, 0.5);
            beta = mvregress(motPreds',trace);
            traceMot = motPreds'*beta;
            corrected = trace' - traceMot';
        end

        function loadROIsDirect (obj, arg1)
            if ~(any(obj.drsave)) || isempty(obj.drsave)
                obj.drsave = obj.dr;
            end
%             [fn dr] = uigetfile([obj.drsave filesep '*_ROIS.mat']);
            load(arg1, 'roiData');

            try
                delete(obj.hROIs);
            catch
            end
            obj.hROIs = [];
            disp(length(roiData))
            for rix = 1:length(roiData)
                %Lucas added this
                switch roiData{rix}.Type
                    case 'images.roi.ellipse'
                        if isempty(roiData{rix}.Center)
                            disp(rix)
                            disp('ellipse was incomplete, creating a dud ellipse to delete later')
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', 0.9, 'FixedAspectRatio', 0, 'Center', [5, 5], 'SemiAxes', [3, 3], 'RotationAngle', 215, 'Label', 'dud');                        
                        else
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', roiData{rix}.AspectRatio, 'FixedAspectRatio', roiData{rix}.FixedAspectRatio, 'Center', roiData{rix}.Center, 'SemiAxes', roiData{rix}.SemiAxes, 'RotationAngle', roiData{rix}.RotationAngle, 'Label', roiData{rix}.Label);
                        end
                    case 'images.roi.circle'
                        if isempty(roiData{rix}.Center)
                            disp(rix)
                            disp('circle was incomplete, creating a dud circle to delete later')
                            obj.hROIs(rix) = images.roi.Circle(obj.hAx,'Center', [5, 5], 'Label', 'dud');
                        else
                            obj.hROIs(rix) = images.roi.Circle(obj.hAx,'Center', roiData{rix}.Center, 'Label', roiData{rix}.Label);
                        end
                    
                    case 'images.roi.polygon'
                        if isempty(roiData{rix}.Position)
                            disp(rix)
                            disp('Polygon was incomplete, creating a dud ellipse to delete later')
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', 0.9, 'FixedAspectRatio', 0, 'Center', [5, 5], 'SemiAxes', [3, 3], 'RotationAngle', 215, 'Label', 'dud');                        
                        else
                            obj.hROIs(rix) = images.roi.Polygon(obj.hAx,'Position', roiData{rix}.Position, 'Label', roiData{rix}.Label);
                        end

                    otherwise
                        error('bad ROI data type');
                end

                fnms = fieldnames(roiData{rix});
                for ix =1:length(fnms)
                    try
                        set(obj.hROIs(rix), fnms{ix}, roiData{rix}.(fnms{ix}));
                    catch
                    end
                end


%                 cm = get(obj.hROIs(rix), 'ContextMenu');
%                 uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(rix)))); %add a 'Set Label' context menu
            end
        end

        function loadROIs (obj, arg1, arg2)
            if ~(any(obj.drsave)) || isempty(obj.drsave)
                obj.drsave = obj.dr;
            end
            [fn dr] = uigetfile([obj.drsave filesep '*_ROIS.mat']);
            load([dr filesep fn], 'roiData');

            try
                delete(obj.hROIs);
            catch
            end
            obj.hROIs = [];
            disp(length(roiData))
            for rix = 1:length(roiData)
                %Lucas added this
                switch roiData{rix}.Type
                    case 'images.roi.ellipse'
                        if isempty(roiData{rix}.Center)
                            disp(rix)
                            disp('ellipse was incomplete, creating a dud ellipse to delete later')
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', 0.9, 'FixedAspectRatio', 0, 'Center', [5, 5], 'SemiAxes', [3, 3], 'RotationAngle', 215, 'Label', 'dud');                        
                        else
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', roiData{rix}.AspectRatio, 'FixedAspectRatio', roiData{rix}.FixedAspectRatio, 'Center', roiData{rix}.Center, 'SemiAxes', roiData{rix}.SemiAxes, 'RotationAngle', roiData{rix}.RotationAngle, 'Label', roiData{rix}.Label);
                        end
                    case 'images.roi.circle'
                        if isempty(roiData{rix}.Center)
                            disp(rix)
                            disp('circle was incomplete, creating a dud circle to delete later')
                            obj.hROIs(rix) = images.roi.Circle(obj.hAx,'Center', [5, 5], 'Label', 'dud');
                        else
                            obj.hROIs(rix) = images.roi.Circle(obj.hAx,'Center', roiData{rix}.Center, 'Label', roiData{rix}.Label);
                        end
                    
                    case 'images.roi.polygon'
                        if isempty(roiData{rix}.Position)
                            disp(rix)
                            disp('Polygon was incomplete, creating a dud ellipse to delete later')
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', 0.9, 'FixedAspectRatio', 0, 'Center', [5, 5], 'SemiAxes', [3, 3], 'RotationAngle', 215, 'Label', 'dud');                        
                        else
                            obj.hROIs(rix) = images.roi.Polygon(obj.hAx,'Position', roiData{rix}.Position, 'Label', roiData{rix}.Label);
                        end

                    otherwise
                        error('bad ROI data type');
                end

                fnms = fieldnames(roiData{rix});
                for ix =1:length(fnms)
                    try
                        set(obj.hROIs(rix), fnms{ix}, roiData{rix}.(fnms{ix}));
                    catch
                    end
                end


                cm = get(obj.hROIs(rix), 'ContextMenu');
                uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(rix)))); %add a 'Set Label' context menu
            end
        end

        function drawPolygon(obj, arg1, arg2)
            L = length(obj.hROIs)+1;
            obj.hROIs(L) = drawpolygon(obj.hAx);
            set(obj.hROIs(L), 'FaceAlpha', 0, 'labelAlpha', 0 , 'labelTextColor', 'y', 'linewidth', 0.5, 'MarkerSize', 1);
            cm = get(obj.hROIs(L), 'ContextMenu');
            uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(L)))); %add a 'Set Label' context menu
            obj.setLabel(obj.hROIs(L));
        end

        function drawCircle(obj, arg1, arg2)
            L = length(obj.hROIs)+1;
            obj.hROIs(L) = drawellipse(obj.hAx);
            set(obj.hROIs(L), 'FaceAlpha', 0, 'labelAlpha', 0 , 'labelTextColor', 'y', 'linewidth', 0.5, 'MarkerSize', 1);
            cm = get(obj.hROIs(L), 'ContextMenu');
            uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(L)))); %add a 'Set Label' context menu
            obj.setLabel(obj.hROIs(L));
        end

        function setLabel(obj, hROI)
            label = inputdlg('Set a label');
            set(hROI, 'Label', label{1});
        end
        function kpf(obj, arg1, arg2)
            switch arg2.Key
                case 'a'
                    @obj.drawPolygon;
                case 'c'
                    @obj.drawCircle;
                case 'period'
                    obj.Clevel = obj.Clevel/1.5;
                    obj.updateDisplay
                case 'comma'
                    obj.Clevel = obj.Clevel.*1.5;
                    obj.updateDisplay;
            end
        end
    end
end