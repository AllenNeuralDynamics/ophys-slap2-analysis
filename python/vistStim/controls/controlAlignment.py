#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import tifffile as tiff
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
slap2ClockPath      = 'Z:/ophys/V1_VisStim/HARP/syncTest/20241008/frameClk_0.csv'
gradientOrderPath   = 'Z:/ophys/V1_VisStim/HARP/syncTest/20241008/orientations_0.csv'
photoDiodeTracePath = 'Z:/ophys/V1_VisStim/HARP/syncTest/20241008/photodiode_0.csv'
expmtFile           = 'Z:/ophys/V1_VisStim/experiments/syncTest_michaelDriftingGratingsSLAP2/'
saveTracesHere      = 'Z:/ophys/V1_VisStim/experiments/controls/'

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

## Extract SLAP2 Signal
'''
# (option A - control) -- from raw tiffs and load into control npy files
dmd1Trials = glob(expmtFile + '*dmd1*.tif')
dmd2Trials = glob(expmtFile + '*dmd2*.tif')
trial=1
for dmd1TrialPath in dmd1Trials[1:]: #ignoring trial 1 (b/c its just spont activity)
    with tiff.TiffFile(dmd1TrialPath) as tif:
        tiffStack = tif.asarray()
    zProfile = np.nanmean(tiffStack, axis=(1,2))
    np.save(saveTracesHere + f'dmd1_{trial}',zProfile[::2])
    trial+=1

trial=1
for dmd2TrialPath in dmd2Trials[1:]: #ignoring trial 1 (b/c its just spont activity)
    with tiff.TiffFile(dmd2TrialPath) as tif:
        tiffStack = tif.asarray()
    zProfile = np.nanmean(tiffStack, axis=(1,2))
    np.save(saveTracesHere+ f'dmd2_{trial}',zProfile[::2])
    trial+=1
'''
dmd1Trials   = glob(saveTracesHere + '*dmd1*')
dmd2Trials   = glob(saveTracesHere + '*dmd2*')
dmd1Trials   = natSort(dmd1Trials)
dmd2Trials   = natSort(dmd2Trials)

## Plot Tiff Traces overlayed on top of Photodiode Signal
plt.close('all')
tix = 0
lastTrial = False
for trialPath in dmd1Trials:
    trialData = np.load(trialPath)
    if tix == 14:
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
