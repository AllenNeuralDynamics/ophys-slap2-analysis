function processSLAP2(dr)

if ~nargin
    dr = uigetdir; %neuron folder where scans are, not project folder
end
aParams = setParams('multiRoiRegSLAP2');
sParams.microscope = "SLAP2";
sParams.drawUserRois = true;
sParams.nParallelWorkers = 8;
sParams = setParams('summarize_NoLoCo', sParams, true);

fullPathToTrialTable = [dr filesep 'trialTable.mat'];
%generate the trial table
if ~exist(fullPathToTrialTable, 'file')
    buildTrialTableSLAP2(dr);
end


%align files
multiRoiRegSLAP2(fullPathToTrialTable,aParams)

%summarize
summarize_NoLoCo(dr, sParams);