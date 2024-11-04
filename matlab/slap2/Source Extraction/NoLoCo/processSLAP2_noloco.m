function processSLAP2_noloco(dr)

if ~nargin
    dr = uigetdir; %neuron folder where scans are, not project folder
end

sParams = setParams('summarizeSLAP2');

%generate the trial table
if ~exist([dr filesep 'trialTable.mat'], 'file')
    buildTrialTable(dr);
end

%align files
multiRoiRegBCI(dr)

%summarize
summarizeSLAP2(dr, sParams);