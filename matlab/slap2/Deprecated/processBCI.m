
function processBCI(dr)

if ~nargin
    dr = uigetdir; %neuron folder where scans are, not project folder
end

%generate the trial table
if ~exist([dr filesep 'trialTable.mat'], 'file')
    buildTrialTable(dr);
end

%align files
multiRoiRegBCI(dr)

%summarize
% summarizeBCI(dr);