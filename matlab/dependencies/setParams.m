function params = setParams(fnName, paramsIn, forceGUI)
%A shared parameter-setting function

switch fnName
    case 'summarize_NoLoCo'
        params.microscope = { '''bergamo''', '''SLAP2'''};          tooltips.scope = 'SLAP2 or bergamo';
        params.includeIntegrationROIs = false; tooltips.includeIntegrationROIs = 'Use integration ROIs for trace extraction?';
        params.sigma_px = 1.33;          tooltips.sigma_px = 'Estimated radius of the PSF (gaussian sigma)';
        params.sparseFac = 0.05;         tooltips.sparseFac = 'sparsity factor for shrinking sources in space, 0-1, higher value makes things sparser';
        params.nmfIter = 5;              tooltips.nmfIter = 'number of iterations of NMF refinement';
        params.dXY = 4;                  tooltips.dXY = 'how large sources can be (radius), pixels';
        params.denoiseWindow_s = 0.25;   tooltips.denoiseWindow_s= 'the timescale on which signals can be smoothed when denoising, seconds';
        params.baselineWindow_Glu_s = 4; tooltips.baselineWindow_Glu_s= 'timescale for calculating F0 in glutamate channel, seconds';
        params.baselineWindow_Ca_s = 4;  tooltips.baselineWindow_Ca_s= 'timescale for calculating F0 in calcium channel, seconds';
        params.activityChannel = 1;      tooltips.activityChannel = 'the channel of the original tiff image that contains the glutamate signal';
        params.tau_s = 0.05;             tooltips.tau_s = 'decay time constant of glutamate signal';
        params.maxSynapseDensity = 0.01; tooltips.maxSynapseDensity = 'maximum synapses per pixel';
        params.nParallelWorkers = 6;     tooltips.nWorkers = 'number of parallel workers';
        params.drawUserRois = false;     tooltips.drawUserRois = 'pop up a GUI to annotate user ROIs?';  
        params.motionThresh = 2.5;       tooltips.motionThresh = 'decrease this to be more stringent on motion correction when censoring frames';
        params.analyzeHz = 200;          tooltips.analyzeHz = 'frame rate used for analysis (SLAP2 only)';
        params.nanThresh = 0.33;         tooltips.nanThresh = 'Max fraction of samples that can be NaN for including a pixel in analysis';
        params.roiHz = 50;               tooltips.roiHz = 'Frame rate for extracting ROI signals, including soma';
        params.discardInitial_s = 0;     tooltips.discardInitial_s = 'time in seconds to remove from analysis at the start of each trial, to accound for warmup';
        params.operator = 'Maria Goeppert Mayer';       tooltips.operator = 'person running the analysis';
        params.makeJSON = false;             tooltips.makeJSON = 'run python script to create processing.json';
    case 'summarize_LoCo'
        params.microscope = { '''SLAP2''' , '''bergamo'''};          tooltips.scope = 'SLAP2 or bergamo';
        params.includeIntegrationROIs = false; tooltips.includeIntegrationROIs = 'Use integration ROIs for trace extraction?';
        params.sigma_px = 1.33;          tooltips.sigma_px = 'Estimated radius of the PSF (gaussian sigma)';
        params.nmfIter = 2;              tooltips.nmfIter = 'number of iterations of NMF refinement';
        params.dXY = 3;                  tooltips.dXY = 'how large sources can be (radius), pixels';
        params.lambda = [];              tooltips.lambda = 'regularizer; roughly the single-photon amplitude. Leave empty to use default/estimate from data.';
        params.denoiseWindow_s = 0.2;   tooltips.denoiseWindow_s= 'the timescale on which signals can be smoothed when denoising, seconds';
        params.baselineWindow_Glu_s = 4; tooltips.baselineWindow_Glu_s= 'timescale for calculating F0 in glutamate channel, seconds';
        params.baselineWindow_Ca_s = 4;  tooltips.baselineWindow_Ca_s= 'timescale for calculating F0 in calcium channel, seconds';
        params.activityChannel = 1;      tooltips.activityChannel = 'the channel of the original tiff image that contains the glutamate signal';
        params.tau_s = 0.03;             tooltips.tau_s = 'decay time constant of glutamate signal';
        params.tau2_s = 0.15;            tooltips.tau2_s = 'decay time constant of 2nd channel signal at synapses (usually spine calcium)';
        params.poissBasedStdIM = 0;      tooltips.poissBasedStdIM = 'use Poisson model to estimate stdIM';
        params.VIF = 1;                  tooltips.VIF = 'variance inflation factor for stdIM estimate (Poiss-based only)';
        params.peakth = 3.5;             tooltips.peakth = 'peak identification threshold (actIM z-score)';
        params.peakFuncOpt = 2;             tooltips.peakFuncOpt = 'peak fitting function (1=gaussian, 2=binned gaussian)';
        params.actImHeteroscedasticNoise = 1;             tooltips.actImHeteroscedasticNoise = 'noise in actIM modeled as heteroscedastic';
        params.peakBufferSize = 0;             tooltips.peakBufferSize = 'number of pixels to mask around each peak';
        params.nParallelWorkers = 12;    tooltips.nWorkers = 'number of parallel workers';
        params.drawUserRois = true;     tooltips.drawUserRois = 'pop up a GUI to annotate user ROIs?';  
        params.motionThresh = 2.5;       tooltips.motionThresh = 'decrease this to be more stringent on motion correction when censoring frames';
        params.analyzeHz = 200;          tooltips.analyzeHz = 'frame rate used for analysis (SLAP2 only)';
        params.nanThresh = 0.33;         tooltips.nanThresh = 'Max fraction of samples that can be NaN for including a pixel in analysis';
        params.discardInitial_s = 0;     tooltips.discardInitial_s = 'time in seconds to remove from analysis at the start of each trial, to accound for warmup';
        params.operator = 'Maria Goeppert Mayer';       tooltips.operator = 'person running the analysis';
        params.makeJSON = false;             tooltips.makeJSON = 'run python script to create processing.json';
    case 'multiRoiRegSLAP2'
        params.alignHz = 80; tooltips.alignHz = 'Frequency for generating downsampled aligned tiffs';
        params.maxshift = 40; tooltips.maxshift = 'Maximum frame offset,in pixels';
        params.clipShift = 5; tooltips.clipShift = 'Maximum allowable shift per frame';
        params.alpha = 0.005; tooltips.alpha = 'exponential decay of template per frame';%exponential time constant for template
        params.nWorkers = 16; tooltips.nWorkers = 'number of parallel workers';
        params.overwriteExisting = false; tooltips.overwriteExisting = 'Realign and overwrite any existing files?';
        params.refStackTemplate = false; tooltips.refStackTemplate = 'Use ref stack as template';
        params.isReVolt = false; tooltips.isReVolt = 'select true for recordings with simultaneous red 1P imaging';
        params.includeIntegrationROIs = false; tooltips.includeIntegrationROIs = 'Use integration ROIs for alignment and TIFF generation?';
        params.operator = 'Maria Goeppert Mayer';       tooltips.operator = 'person running the analysis';
    case 'integrationRegistration'
        params.alignHz = 80; tooltips.alignHz = 'Frequency for generating downsampled aligned tiffs';
        params.maxshiftXY = 25; tooltips.maxshift = 'Maximum frame offset,in pixels';
        params.maxshiftZ = 10; tooltips.maxshift = 'Maximum frame offset,in pixels';
        params.clipShift = 5; tooltips.clipShift = 'Maximum allowable shift per frame';
        params.motionMetric = {'''poisson''','''correlation'''}; tooltips.motionMetric = 'Metric for selecting best motion shift';
        params.robust = false; tooltips.robust = 'Use robust likelihood?';
        params.efficientTiffSave = false; tooltips.efficientTiffSave = 'Save Tiffs locally first then transfer?';
        params.tempFileDir = 'C:\temp'; tooltips.tempFileDir = 'Directory for temp files';
        % params.alpha = 0.005; tooltips.alpha = 'exponential decay of template per frame';%exponential time constant for template
        params.nWorkers = 16; tooltips.nWorkers = 'number of parallel workers';
        params.overwriteExisting = false; tooltips.overwriteExisting = 'Realign and overwrite any existing files?';
        params.integrationOnly = false; tooltips.integrationOnly = 'Align only on integration superpixels';
        params.saveTiffs = true; tooltips.saveTiffs = 'Save aligned tiff movies';
        params.operator = 'Maria Goeppert Mayer';       tooltips.operator = 'person running the analysis';
    case 'stripRegBergamo'
        params.maxshift = 50; tooltips.maxshift = 'Maximum frame offset,in pixels';
        params.clipShift = 10; tooltips.clipShift = 'Maximum allowable shift per frame';
        params.nWorkers = 4; tooltips.nWorkers = 'number of parallel workers';
        params.overwriteExisting = false; tooltips.overwriteExisting = 'Realign and overwrite any existing files?';
        params.removeLines = 4; tooltips.removeLines = 'remove this many flyback lines from the top of each image';
        params.ds_time = 3; tooltips.ds_time = 'movies are downsampled (2^ds_time)x in time for alignment';
        params.frameRate = 0; tooltips.frameRate = 'imaging frame rate; if 0, calculated from metadata or set as default';
        params.denoise20Hz = false;
        params.saveTif = true; tooltips.saveTif = 'whether to save registered movie as .tif or .h5';
    case {'extractDendrites', 'extractDendrites_new'}
        % Epoch-aware SLAP2 dendritic-voltage extraction. extractDendrites_new
        % remains accepted for backward compatibility, but extractDendrites is
        % now the canonical function name.
        params.manualROIs = false;
        tooltips.manualROIs = 'Draw ROIs manually? If false, use SLAP2 integration ROIs.';

        params.chIdx = 1;
        tooltips.chIdx = 'Channel index passed to the SLAP2 Trace object.';

        params.zIdx = 1;
        tooltips.zIdx = 'Z-plane index passed to the SLAP2 Trace object. Usually 1 for single-plane/DMD voltage extraction.';

        params.windowWidth_lines = 16;
        tooltips.windowWidth_lines = 'Window width, in SLAP2 lines, passed to Trace.process. Matches the original extractDendrites default.';

        params.expectedWindowWidth_lines = 5000;
        tooltips.expectedWindowWidth_lines = 'Expected/baseline window width, in SLAP2 lines, passed to Trace.process. Matches the original extractDendrites default.';

        params.outputMode = 'trial';
        tooltips.outputMode = 'Trace organization: trial, continuous, or both. Trial is recommended for multi-trial and multi-epoch voltage analysis.';

        params.storageMode = 'h5';
        tooltips.storageMode = 'Storage backend for large trace arrays. h5 is recommended for large recordings.';

        params.precision = 'single';
        tooltips.precision = 'Numeric precision for saved traces. single reduces memory and disk use relative to double.';

        params.outputDir = '';
        tooltips.outputDir = 'Custom output directory. Leave empty to save beside trialTable.mat or in the configured subfolder.';

        params.useOutputSubfolder = true;
        tooltips.useOutputSubfolder = 'Save outputs into a dedicated subfolder beside trialTable.mat. Recommended.';

        params.outputSubfolderName = 'dendriticVoltageExtraction';
        tooltips.outputSubfolderName = 'Name of the output subfolder created beside trialTable.mat.';

        params.timestampOutputSubfolder = false;
        tooltips.timestampOutputSubfolder = 'Create a timestamped output subfolder for each extraction run.';

        params.useParallel = true;
        tooltips.useParallel = 'Use MATLAB parallel workers for ROI extraction. Disable for debugging or lowest-memory runs.';

        params.numWorkers = 4;
        tooltips.numWorkers = 'Requested number of MATLAB parallel workers. 2-4 is usually safer than 16-24 for this I/O-heavy workflow.';

        params.restartPool = true;
        tooltips.restartPool = 'Restart the parallel pool if its worker count does not match numWorkers.';

        params.maxConcurrentROIs = 2;
        tooltips.maxConcurrentROIs = 'Maximum number of ROIs extracted at the same time per DMD. Main memory-throttling parameter.';

        params.makePlots = false;
        tooltips.makePlots = 'Generate ROI/reference preview plots during extraction. Disable for batch processing.';

        params.saveRefImages = true;
        tooltips.saveRefImages = 'Save selected reference images in the lightweight summary .mat file.';

        params.saveSummaryAfterEachBatch = true;
        tooltips.saveSummaryAfterEachBatch = 'Checkpoint the summary .mat file after each ROI batch.';

        params.resume = false;
        tooltips.resume = 'Reserved for resumable extraction. Currently prevents accidental reuse of old HDF5 outputs.';

        params.stopOnError = true;
        tooltips.stopOnError = 'Stop immediately when an ROI fails. Set false to log errors and continue.';

        params.h5ChunkLines = 100000;
        tooltips.h5ChunkLines = 'HDF5 chunk size along the line/time dimension. 50000-200000 is a reasonable range.';

        params.h5Deflate = 0;
        tooltips.h5Deflate = 'HDF5 compression level. 0 is fastest; 1-9 saves space but slows writing.';
    case 'summarize_Voltage'
        params.tau = 0.001;            tooltips.tau_s = 'decay time constant of voltage signal';
        params.analyzeHz = 1000; %frame rate used for analysis
        params.discardInitial_s = 0.1; %discard the first short period of each trial as the beam stabilization locks on and the imaging system warms up
        params.denoiseWindow_s = 0.25; %number of samples to average together for denoising
        params.baselineWindow_s = 1; %timescale for calculating F0 in glutamate channel, seconds
        params.exptType = {'"V1 Gratings"', '"BCI"', '"other"'};  tooltips.exptType ='Experiment type';
        params.motionThresh = 2.5;       tooltips.motionThresh = 'decrease thresh to be more stringent on motion correction when censoring frames';   
        params.manualRois = false;      tooltips.manualRois = 'Redraw ROIs manually? If false, will use the ROIs drawn at imaging time';
        params.integrationOnly = true;  tooltips.integrationOnly = 'Only process integration ROIs?';
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