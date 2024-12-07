function params = setParams(fnName, paramsIn, forceGUI)
%A shared parameter-setting function

switch fnName

    case 'summarizeBergamo_NoLoCo'
        params.sigma_px = 1.33;          tooltips.sigma_px = 'Estimated radius of the PSF (gaussian sigma)';
        params.sparseFac = 0.05;         tooltips.sparseFac = 'sparsity factor for shrinking sources in space, 0-1, higher value makes things sparser';
        params.nmfIter = 5;              tooltips.nmfIter = 'number of iterations of NMF refinement';
        params.dXY = 4;                  tooltips.dXY = 'how large sources can be (radius), pixels';
        params.denoiseWindow_s = 0.25;   tooltips.denoiseWindow_s= 'the timescale on which signals can be smoothed when denoising, seconds';
        params.baselineWindow_Glu_s = 4; tooltips.baselineWindow_Glu_s= 'timescale for calculating F0 in glutamate channel, seconds';
        params.baselineWindow_Ca_s = 4;  tooltips.baselineWindow_Ca_s= 'timescale for calculating F0 in calcium channel, seconds';
        params.activityChannel = 1;      tooltips.activityChannel = 'which channel index in the tiff image contains the glutamate signal (ch1 is the first channel recorded, etc regardless of the scanimage digitizer channel)';
        params.tau_s = 0.05;             tooltips.tau_s = 'decay time constant of glutamate signal';
        params.maxSynapseDensity = 0.01; tooltips.maxSynapseDensity = 'maximum synapses per pixel';
        params.motionThresh = 2;         tooltips.motionThresh = 'decrease thresh to be more stringent on motion correction when censoring frames';
        params.nanThresh = 0.25;        tooltips.nanThresh = 'Max fraction of samples that can be NaN for including a pixel in analysis';
    case 'summarize_NoLoCo'
        params.microscope = { '"bergamo"', '"SLAP2"'};          tooltips.scope = 'SLAP2 or bergamo';
        params.sigma_px = 1.33;          tooltips.sigma_px = 'Estimated radius of the PSF (gaussian sigma)';
        params.sparseFac = 0.05;         tooltips.sparseFac = 'sparsity factor for shrinking sources in space, 0-1, higher value makes things sparser';
        params.nmfIter = 5;              tooltips.nmfIter = 'number of iterations of NMF refinement';
        params.dXY = 4;                  tooltips.dXY = 'how large sources can be (radius), pixels';
        params.denoiseWindow_s = 0.25;   tooltips.denoiseWindow_s= 'the timescale on which signals can be smoothed when denoising, seconds';
        params.baselineWindow_Glu_s = 4; tooltips.baselineWindow_Glu_s= 'timescale for calculating F0 in glutamate channel, seconds';
        params.baselineWindow_Ca_s = 4;  tooltips.baselineWindow_Ca_s= 'timescale for calculating F0 in calcium channel, seconds';
        params.activityChannel = 2;      tooltips.activityChannel = 'the channel of the original tiff image that contains the glutamate signal';
        params.tau_s = 0.05;             tooltips.tau_s = 'decay time constant of glutamate signal';
        params.maxSynapseDensity = 0.01; tooltips.maxSynapseDensity = 'maximum synapses per pixel';
        params.nParallelWorkers = 8;     tooltips.nWorkers = 'number of parallel workers';
        params.drawUserRois = false;      tooltips.drawUserRois = 'pop up a GUI to annotate user ROIs?';  
        params.motionThresh = 2;         tooltips.motionThresh = 'decrease this to be more stringent on motion correction when censoring frames';
        params.analyzeHz = 200;          tooltips.analyzeHz = 'frame rate used for analysis (SLAP2 only)';
        params.nanThresh = 0.25;        tooltips.nanThresh = 'Max fraction of samples that can be NaN for including a pixel in analysis';
    case 'summarizeSLAP2'
        params.tau_s = 0.05;            tooltips.tau_s = 'decay time constant of glutamate signal';
        params.analyzeHz = 200; %frame rate used for analysis
        params.discardInitial_s = 0.1; %discard the first short period of each trial as the beam stabilization locks on and the imaging system warms up
        params.sigma_px = 1.5; tooltips.sigma_px = 'Estimated radius of the PSF (gaussian sigma)';
        params.denoiseWindow_s = 0.25; %number of samples to average together for denoising
        params.baselineWindow_Glu_s = 4; %timescale for calculating F0 in glutamate channel, seconds
        params.baselineWindow_Ca_s = 4; %timescale for calculating F0 in calcium channel, seconds
        params.sparseFac = 0.05; %sparsity factor for shrinking sources in space, 0-1, higher value makes things sparser
        params.nmfIter = 4; %number of iterations of NMF refinement
        params.dXY = 5; %how large sources can be (radius), pixels
        params.exptType = {'"V1 Gratings"', '"BCI"', '"other"'};  tooltips.exptType ='Experiment type';
        params.maxSynapseDensity = 0.01; tooltips.maxSynapseDensity = 'maximum synapses per pixel';
        params.motionThresh = 2.5;       tooltips.motionThresh = 'decrease thresh to be more stringent on motion correction when censoring frames';
    case 'multiRoiRegSLAP2'
        params.alignHz = 80; tooltips.alignHz = 'Frequency for generating downsampled aligned tiffs';
        params.maxshift = 50; tooltips.maxshift = 'Maximum frame offset,in pixels';
        params.clipShift = 5; tooltips.clipShift = 'Maximum allowable shift per frame';
        params.alpha = 0.005; tooltips.alpha = 'exponential decay of template per frame';%exponential time constant for template
        params.nWorkers = 16; tooltips.nWorkers = 'number of parallel workers';
        params.overwriteExisting = false; tooltips.overwriteExisting = 'Realign and overwrite any existing files?';
        params.refStackTemplate = false; tooltips.refStackTemplate = 'Use ref stack as template';
    case 'integrationRegistration'
        params.alignHz = 80; tooltips.alignHz = 'Frequency for generating downsampled aligned tiffs';
        params.maxshiftXY = 25; tooltips.maxshift = 'Maximum frame offset,in pixels';
        params.maxshiftZ = 10; tooltips.maxshift = 'Maximum frame offset,in pixels';
        params.clipShift = 5; tooltips.clipShift = 'Maximum allowable shift per frame';
        params.robust = false; tooltips.robust = 'Use robust likelihood?';
        params.efficientTiffSave = false; tooltips.efficientTiffSave = 'Save Tiffs locally first then transfer?';
        params.tempFileDir = 'C:\temp'; tooltips.tempFileDir = 'Directory for temp files';
        % params.alpha = 0.005; tooltips.alpha = 'exponential decay of template per frame';%exponential time constant for template
        params.nWorkers = 16; tooltips.nWorkers = 'number of parallel workers';
        params.overwriteExisting = false; tooltips.overwriteExisting = 'Realign and overwrite any existing files?';
    case 'stripRegBergamo'
        params.maxshift = 50; tooltips.maxshift = 'Maximum frame offset,in pixels';
        params.clipShift = 10; tooltips.clipShift = 'Maximum allowable shift per frame';
        params.nWorkers = 4; tooltips.nWorkers = 'number of parallel workers';
        params.overwriteExisting = false; tooltips.overwriteExisting = 'Realign and overwrite any existing files?';
        params.removeLines = 4; tooltips.removeLines = 'remove this many flyback lines from the top of each image';
        params.ds_time = 3; tooltips.ds_time = 'movies are downsampled (2^ds_time)x in time for alignment';
        params.denoise20Hz = false;
    otherwise
        error('Unknown function name passed to setParams.m')
end             

if nargin>1 %if the user specified parameters, add in the user parameters, use defaults for remaining, NO GUI
    for field = fieldnames(paramsIn)'
         params.(field{1}) = paramsIn.(field{1});
    end
    if nargin<3 || ~forceGUI
        return
    end
end

%get parameters from user
% prefPath = [fileparts(which(fnName)) filesep fnName '_prefs.mat'];
% % if exist(prefPath, 'file')
% %     paramsIn = load(prefPath, 'paramsIn');
% %     % if (length(fieldnames(paramsIn)) == length(fieldnames(params))) && all(strcmp(fieldnames(paramsIn), fieldnames(params)))
% %     %     params = paramsIn;
% %     % end
% % end
 paramsIn = optionsGUI(params, tooltips, fnName);
 params = paramsIn;
% save(prefPath, 'paramsIn');

end