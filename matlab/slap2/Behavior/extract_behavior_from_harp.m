function behaviorData = extract_behavior_from_harp(dr)
    %{
        given a slap2 data dir, find the appropriate behavior dir
        SUBJECT TO CHANGE DEPENDING ON SCI-COM CONCENSUS FOR DATA ORGANIZATION
    
    Inputs:
        dr (str) : path to SLAP2 Directory
    Outputs:
        behaviorData (cell structure array) : cell array of trial structures
        containing extracted Harp data from each trial
        
       *****NOTE*****************************************************
       *cell structure arr is organized by epoch, should be         *
       *indexable based on initial epoch selection in processBCI    *
       **************************************************************
    %}
    
    mousePath = [dr(1:(strfind(dr, 'SLAP2')-1)) 'Behavior\BCI\'];
    mouseBehaviorDr = dir(mousePath);
    temp = split(dr, '\'); sessionDate_string = temp{end-1};
    temp2 = split(sessionDate_string, '_'); dateID = erase(temp2{3}, '-');
    epochIndexing = cell2mat(arrayfun(@(i) contains(mouseBehaviorDr(i).name, dateID), 1:numel(mouseBehaviorDr(:)), 'UniformOutput', false));
    
    mousePaths = ls(mousePath);
    epochsOfInterest = mousePaths(epochIndexing,:);
    epochPaths = arrayfun(@(i) [mousePath epochsOfInterest(i,:)], 1:size(epochsOfInterest,1), 'UniformOutput',false);
    epochTypes = contains(epochPaths, 'BCIon');
    moveOnly = epochPaths{~epochTypes};
    BCIon = epochPaths{epochTypes}; BCIon = BCIon(find(~isspace(BCIon)));
    
    registers = {
        '\Behavior.harp\Register__PwmFrequencyDO1.bin', ...
        '\Behavior.harp\Register__OutputSet.bin',...
        '\Behavior.harp\Register__DigitalInputState.bin',...
        '\Behavior.harp\Register__AnalogData.bin',...
        '\Operation\SpoutPosition.csv',...
        '\LoadCells.harp\Register__LoadCellData.bin'
        }; %registers subject to change based on updates
    
    %HARP PINOUT
    startPin = 0x1000; stopPin = 0x1; %Subject to change based on bonsai updates
    goCue = 500; rewardCue = 100;
    licks = 2;
    
    %MoveOnlys
    moveOnly_audio = read_harp_bin([moveOnly registers{1}]);
    moveOnly_goCueT = moveOnly_audio.Seconds(moveOnly_audio.Values==goCue);
    moveOnly_rewCueT = moveOnly_audio.Seconds(moveOnly_audio.Values==rewardCue);
    
    moveOnly_handshake = read_harp_bin([moveOnly registers{2}]);
    moveOnly_startTriggers = moveOnly_handshake.Seconds(moveOnly_handshake.Values == startPin);
    moveOnly_stopTriggers = moveOnly_handshake.Seconds(moveOnly_handshake.Values == stopPin);
    moveOnly_stopTriggers = moveOnly_stopTriggers(2:end);%wierd instance of a stop trigger being sent before the first start trigger
    moveOnly_trialTimes = moveOnly_stopTriggers - moveOnly_startTriggers;
    
    
    moveOnly_licks = read_harp_bin([moveOnly registers{3}]);
    moveOnly_lickT = moveOnly_licks.Seconds(moveOnly_licks.Values==licks);
    
    moveOnly_adc = read_harp_bin([moveOnly registers{4}]); moveOnly_adc.Values = moveOnly_adc.Values(:,1);
    moveOnly_motor = readtable([moveOnly registers{5}]);
    moveOnly_loadCell = read_harp_bin([moveOnly registers{6}]);  
    
    moveOnly_array = trialize_data(moveOnly_startTriggers, moveOnly_stopTriggers, moveOnly_loadCell, moveOnly_adc, moveOnly_motor, moveOnly_lickT, moveOnly_goCueT, moveOnly_rewCueT);
    
    %BCIs
    BCI_audio = read_harp_bin([BCIon registers{1}]);
    BCI_goCueT = BCI_audio.Seconds(BCI_audio.Values==goCue);
    BCI_rewCueT = BCI_audio.Seconds(BCI_audio.Values==rewardCue);
    
    BCI_handshake = read_harp_bin([BCIon registers{2}]);
    BCI_startTriggers = BCI_handshake.Seconds(BCI_handshake.Values == startPin);
    BCI_stopTriggers = BCI_handshake.Seconds(BCI_handshake.Values == stopPin);
    BCI_stopTriggers = BCI_stopTriggers(2:end);%wierd instance of a stop trigger being sent before the first start trigger
    BCI_trialTimes = BCI_stopTriggers - BCI_startTriggers;
    
    BCI_licks = read_harp_bin([BCIon registers{3}]);
    BCI_lickT = BCI_licks.Seconds(BCI_licks.Values==licks);
    
    BCI_adc = read_harp_bin([BCIon registers{4}]); moveOnly_adc.Values = moveOnly_adc.Values(:,1);
    BCI_motor = readtable([BCIon registers{5}]);
    BCI_loadCell = read_harp_bin([BCIon registers{6}]);
    
    BCI_array = trialize_data(BCI_startTriggers, BCI_stopTriggers, BCI_loadCell, BCI_adc, BCI_motor, BCI_lickT, BCI_goCueT, BCI_rewCueT);
    
    
    behaviorData = vertcat(moveOnly_array, BCI_array);  


end