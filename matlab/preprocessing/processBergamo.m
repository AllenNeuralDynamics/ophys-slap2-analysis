function processBergamo(dr)

if ~nargin
    dr = uigetdir; %neuron folder where scans are, not project folder
end
aParams = setParams('stripRegBergamo');

sParams.microscope = "bergamo";
sParams.drawUserRois = false;
sParams.nParallelWorkers = 3;
sParams = setParams('summarize_NoLoCo', sParams, true);

fullPathToTrialTable = [dr filesep 'trialTable.mat'];
%generate the trial table
if ~exist(fullPathToTrialTable, 'file')
    buildTrialTableBergamo(dr);
end

%align files
stripRegBergamo(dr, 'trialTable.mat',aParams)

%summarize
summarize_NoLoCo(dr, sParams);