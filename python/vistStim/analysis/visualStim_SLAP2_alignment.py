#%%
import mat73, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
from glob import glob
bks = matplotlib.rcsetup.interactive_bk
plt.switch_backend('WebAgg')

## Prelim Functions
def natSort(sList):
    numbersOnly = []
    for s in sList:
        latter = s.split('_')[2]
        number = int(latter.split('.')[0])
        numbersOnly.append(number)
    a = np.argsort(numbersOnly)
    sList = np.array(sList)
    return sList[a[:]]


## Set Paths to Data
slap2ClockPath      = 'Z:/ophys/V1_VisStim/HARP/760268/20241007/frameClk_0.csv'
gradientOrderPath   = 'Z:/ophys/V1_VisStim/HARP/760268/20241007/orientations_0.csv'
photoDiodeTracePath = 'Z:/ophys/V1_VisStim/HARP/760268/20241007/photodiode_0.csv'
expmtFile           = 'Z:/ophys/V1_VisStim/experiments/slap2_760268_2024-10-07_13-35-05/Neuron1/ExperimentSummary/Summary-241008-085607.mat'
saveTracesHere      = 'Z:/ophys/V1_VisStim/experiments/slap2_760268_2024-10-07_13-35-05/Neuron1/ExperimentSummary/traces/'
# if ~os.path.exists(saveTracesHere):
#     os.mkdir(saveTracesHere)
## Read In Data From HARP/Bonsai
photoDiode   = pd.read_csv(photoDiodeTracePath)
slap2Clock   = pd.read_csv(slap2ClockPath)
slap2Clock   = slap2Clock.rename(columns=
                  {slap2Clock.columns[0]:'bools', slap2Clock.columns[1]:'ticks'})
ledTrace     = photoDiode['Item2'].values
ledTimes     = photoDiode['Item1'].values - photoDiode['Item1'].values[0]

## Find offset between aquisition and photodiode starts - use to generate trial starts and stops
offset       =  slap2Clock['ticks'].values[0] - photoDiode['Item1'].values[0]
time         =  (slap2Clock['ticks'].values - slap2Clock['ticks'].values[0])
trialEdges   =  np.diff(time)>0.06
timeMinusOne =  time[:-1]
timePlusOne  =  time[1:]
trialEndings =  timeMinusOne[trialEdges] + offset 
trialStarts  =  timePlusOne[trialEdges]  + offset

##############################################
## Extract SLAP2 Signal From ExpmtSummary
# tempMat = mat73.loadmat(expmtFile, use_attrdict=True)
# ## Somatic Calcium
# for trial in range(len(tempMat['exptSummary']['E'])-1):
#     trialData               = tempMat['exptSummary']['E'][trial+1]
#     dmd1                    = trialData[0]
#     dmd2                    = trialData[1]
#     calciumTracePerTrial    = dmd1['ROIs']['F'][:,1] #red channel
    
#     np.save(saveTracesHere+f'soma_ca_{trial}', np.array(calciumTracePerTrial)) #temporarily save as npy files for easy work later
##############################################

somaTrials =  glob(saveTracesHere + '*soma*')
## Plot Tiff Traces overlayed on top of Photodiode Signal
plt.close('all')
tix = 0
lastTrial = False
for trialPath in somaTrials:
    trialData = np.load(trialPath)
    if tix == 14: #hardcoding last trial...
        lastTrial=True
    if lastTrial:
        startTime       = trialEndings[tix]
        greaterThanT0   = np.where(ledTimes>startTime, True, False)
        ledIDX          = ledTrace[np.array(greaterThanT0 )]
    else:
        startTime       = trialStarts[tix]
        stopTime        = trialEndings[tix+1]
        greaterThanT0   = np.where(ledTimes>startTime, True, False)
        lessThanTf      = np.where(ledTimes<stopTime, True, False)
        ledIDX          = ledTrace[np.array(greaterThanT0 & lessThanTf)]
    harpXaxis   = np.linspace(startTime, stopTime, num=len(ledIDX))
    slap2Xaxis  = np.linspace(startTime, stopTime, num=len(trialData))
    tix += 1
    plt.figure(tix)
    plt.plot(harpXaxis, ledIDX/np.nanmax(ledIDX), color='blue')
    plt.plot(slap2Xaxis, trialData/np.nanmax(trialData),color='orange')
plt.show()


# %%
