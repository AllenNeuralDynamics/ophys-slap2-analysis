function trializedStructure = trialize_data(startTriggers, stopTriggers, loadCellData, adc, motor, licks, gos, rewards)
trializedStructure = struct(); 
trials = cell(size(stopTriggers,1),1);

%going based off final trial because there will always be an end but
%sometimes starts might not have ends
for t = 1:size(stopTriggers,1)
    trialData = struct();
    startTrigger = startTriggers(t);
    stopTrigger = stopTriggers(t);
    

    loadCellTrialTimes = loadCellData.Seconds((loadCellData.Seconds > startTrigger) & (loadCellData.Seconds < stopTrigger));
    loadCellTrialValues = loadCellData.Values((loadCellData.Seconds > startTrigger) & (loadCellData.Seconds < stopTrigger));

    adcTrialTimes = adc.Seconds( (adc.Seconds > startTrigger) & (adc.Seconds < stopTrigger)  );
    adcTrialValues = adc.Values( (adc.Seconds > startTrigger) & (adc.Seconds < stopTrigger)  );

    motorTrialTimes = motor.Seconds(  (motor.Seconds > startTrigger) & (motor.Seconds < stopTrigger) );
    motorTrialValues = motor.Value(  (motor.Seconds > startTrigger) & (motor.Seconds < stopTrigger) );

    lickTimes = licks( (licks > startTrigger) & (licks < stopTrigger));
    goCues = gos( (gos > startTrigger) & (gos < stopTrigger));
    rewardCues = rewards( (rewards > startTrigger) & (rewards < stopTrigger));


    trialData.loadCell.Seconds = loadCellTrialTimes;
    trialData.loadCell.Values = loadCellTrialValues;
    trialData.adc.Seconds = adcTrialTimes;
    trialData.adc.Values = adcTrialValues;
    trialData.Motor.Seconds = motorTrialTimes;
    trialData.Motor.Values = motorTrialValues;
    trialData.licks.Seconds = lickTimes;
    trialData.goCue.Seconds = goCues;
    trialData.rewardCue = rewardCues;
    trialData.rxnTime = stopTrigger - startTrigger;

    trials{t} = trialData;
    
end

trializedStructure = trials;

end
