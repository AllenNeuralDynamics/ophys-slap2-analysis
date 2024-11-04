function params = setParams(fnName, paramsIn)
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
        params.activityChannel = 2;      tooltips.activityChannel = 'the channel of the original tiff image that contains the glutamate signal';
        params.tau_s = 0.027;            tooltips.tau_s = 'decay time constant of glutamate signal';
    case 'summarize_NoLoCo'
        params.scope = 'SLAP2';          tooltips.scope = 'SLAP2 or bergamo';
        params.sigma_px = 1.33;          tooltips.sigma_px = 'Estimated radius of the PSF (gaussian sigma)';
        params.sparseFac = 0.05;         tooltips.sparseFac = 'sparsity factor for shrinking sources in space, 0-1, higher value makes things sparser';
        params.nmfIter = 5;              tooltips.nmfIter = 'number of iterations of NMF refinement';
        params.dXY = 4;                  tooltips.dXY = 'how large sources can be (radius), pixels';
        params.denoiseWindow_s = 0.25;   tooltips.denoiseWindow_s= 'the timescale on which signals can be smoothed when denoising, seconds';
        params.baselineWindow_Glu_s = 4; tooltips.baselineWindow_Glu_s= 'timescale for calculating F0 in glutamate channel, seconds';
        params.baselineWindow_Ca_s = 4;  tooltips.baselineWindow_Ca_s= 'timescale for calculating F0 in calcium channel, seconds';
        params.activityChannel = 2;      tooltips.activityChannel = 'the channel of the original tiff image that contains the glutamate signal';
        params.tau_s = 0.027;            tooltips.tau_s = 'decay time constant of glutamate signal';
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
    otherwise
        error('Unknown function name passed to setParams.m')
end             

if nargin>1 %if the user specified parameters, add in the user parameters, use defaults for remaining, NO GUI
    for field = fieldnames(paramsIn)'
         params.(field{1}) = paramsIn.(field{1});
    end
    return
end

%get parameters from user
prefPath = [fileparts(which(fnName)) filesep fnName '_prefs.mat'];
% if exist(prefPath, 'file')
%     paramsIn = load(prefPath, 'paramsIn');
%     % if (length(fieldnames(paramsIn)) == length(fieldnames(params))) && all(strcmp(fieldnames(paramsIn), fieldnames(params)))
%     %     params = paramsIn;
%     % end
% end
paramsIn = optionsGUI(params, tooltips, fnName);
params = paramsIn;
save(prefPath, 'paramsIn');

end