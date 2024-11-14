## Glutamate data preprocessing

The top-level functions that process glutamate imaging data are `processSLAP2.m` and `processBergamo.m` depending on which microscope the data were acquired on. These functions generate output in a shared format. They can be called with no arguments, in which case a GUI will gather user input. The pipeline consists of three steps, which the top-level function calls in sequence:
1. Defining the trial structure
    * `buildTrialTableSLAP2` or `builtTrialTableBergamo.m`
2. Registering the images within each recording:  
    * ``multiROIRegSLAP2`` or ``stripRegBergamo.m``
3. Extracting signals and generating a "Summary" file for the experiment:
    * `summarize_NoLoCo.m`
  
If processing aborts, you can call the top-level function again and it will skip over any completed steps (unless you set an option flag to repeat them).

When you call these functions without parameters, they pop up a GUI with default parameters that you can change. The default parameters for each function are set in the function `setParams.m`. There are tooltips that explain each parameter when you hover over it.

### Step 1. Defining the trial structure
This step asks the user to identify the files that are part of the experiment, and builds a **trial table** that keeps track of files generated during processing. The experiment is assumed to consist of one or more **epochs**, each consisting of one or more **trials**. All recordings in an experiment are assumed to image the same field of view, and must be in the same folder. The recordings will be aligned to each other and share a common set of extracted synapses; if you image different fields of view in a single session, these should be split into separate folders and processed separately. 

For SLAP2 recordings, the corresponding SLAP2 reference stack should exist in the same folder as the activity recording, or in a subfolder. Do not include multiple reference stacks in the folder; this will cause an error.

When you run this function, a list of files will be shown, and you will be asked to select files belonging to each epoch, in sequence. Files not assigned to an epoch will not be analyzed.

**Bergamo recordings:** If all of your files belong to a single epoch you can automatically generate the trial table by directly calling `stripRegBergamo.m`. This simplifies batch analysis without GUI input.

### Step 2. Registering Images
Image registration aligns the images within each trial to each other and generates a downsampled, aligned copy of each recording as a `.tif` file. It also generates an associated `..._ALIGNMENTDATA.mat` file for each `.tif`, which contains data related to the alignment process and is used by later steps.  

### Step 3. Extracting signals
We represent signals with a matrix factorization model (_i.e._ as a set of spatial patterns that each vary in brightness over time). This approach is shared with common approaches to analyze somatic calcium imaging data (_e.g._ using CaImAn or Suite2P). Our approach has been tailored specifically for synaptic imaging with the following assumptions:
* Sources are small and round (with characteristic size `sigma_px`)
* Sources flash with positive-going transients that rise quickly and have exponential decay time constant `tau_s`.
* Sources frequently overlap with neighboring sources, requiring a superresolution approach to separate them (_e.g._ SOFI or localization microscopy).

Usin the downsampled recordings from step 2, for each trial, we filter the movie in space and time to emphasize sources that match the priors above, subtract out any slow fluctuations, and  compute a high-order statistic (similar to SOFI) to produce an image that sharpens and enhances fluctuating point sources. This image `actIM`, is the activity image for each trial. This computation is done by the function `localizeSources.m` 

We then align all trials in the experiment, averaging together the activity images to create a single activity image for the experiment. Putative synapses are identified as local maxima in this image.

Still using the downsampled recordings, we then form a [pixels x time] matrix from the pixels surrounding the putative synapses across all recordings, and factorize this using constrained nonnegative matrix factorization while enforcing our spatial and temporal priors. During this process, we merge sources that overlap and do not show sufficiently distinct activity (using the Akaike Information Criterion) and discard sources that do not explain sufficient variance in the recordings. 

Finally, we use the resulting footprints of these sources to extract activity from the original, high-framerate trial recordings. The resulting data are saved in an `experimentSummary` file (the filename has the time of processing appended, so you can easily generate and compare versions for different parameter settings).

## The experimentSummary file
The experiment summary file contains a single matlab structure `exptSummary`, with the following fields:
* `E`: A cell array (#FOVs x # Trials) of structures containing the traces and other information extracted from each trial
* `meanIM`: the average image of the sample across all trials, to which each trial is aligned
* `actIM`: the activity image used to identify synatic sites (an average of the aligned per-trial activity images)
* `params`: the user-settable parameters used to process the data

### Information in exptSummary.E{}
Each element of exptSummary.E{} contains data from a particular field of view on a particular trial. It contains the following fields:

* dF ('delta F') : Baseline-subtracted fluorescence
* F0: An estimate of the fluorescence baseline
* dFF ('delta F over F') : Baseline-subtracted fluorescence divided by baseline, i.e. fractional change in fluorescence.
* ROIs : the activity of user-defined ROIs, such as the soma
```
Tip:  
dF vs dFF
When performing analyses of glutamate imaging, the more useful variable is dF, which has units proportional to photons. 

For example, a good estimate of the total measured glutamate input to the cell is the sum of dF across synapses. It is inappropriate to sum dFF in the same way, because very dim synapses can have large dFF values and background fluorescence (e.g. from the dendritic shaft) strongly affects fractional changes at some synapses compared to others. This is different than population Calcium imaging, which tends to work with dF/F. 

In glutamate imaging, dFF is most useful for compensating for bleaching over time for a single synapse, for example to see if its activity increases vs decreases over the course of an experiment. However even for this some confounds remain due to brightness-dependent bias in F0 estimation. 
```
For SLAP2, there are two simultaneously-imaged fields of view, and the traces in the two trials should have the same or very similar numbers of timepoints. For analysis, the two FOVs can be merged- the data are interpolated onto the same timebase.

dF and dFF are themselves structures, with fields corresponding to different versions of the signal with different temporal priors:  
* `ls`: Least Squares solve for the activity of each ROI. No temporal prior or nonnegativity is imposed.
* `matchFilt`: Least squares solve, after applying a matched filter in time. This is the optimal linear filter for detecting isolated release events, and can be thresholded to perform template matching.
* `nonneg`: Generated from matchFilt by deconvolving out the matched filter after solving, with a mild nonnegativity/lower bound constraint. This has only a weak temporal and weak nonnegativity prior, and is similar to ls.  
* `spikes`: Generated from matchFilt by deconvolving out both the matched filter and the expected temporal decay of the indicator, generating sharp spikes at the onsets of events. This removes the effect of indicator kinetics and is useful as an estimate of event times, e.g. to compute sharper crosscorrelation. A weak nonnegativity/lower bound constraint is imposed.
* `denoised`: Generated from matchFilt by re-convolving spikes with the indicator response. This imposes a strong temporal prior and a weak nonnegativity/lower bound constraint. 

```
Tip:
Parallel Processing and Memory usage.
Processing is done in parallel over trials, and the amount of memory needed is proportional to the number of parallel workers and the size of each trial file. For our workstations (128 Gb RAM), the pipeline works best when each downsampled aligned trial is <4 Gb. This allows the computer to use many parallel workers (~15 for 4Gb files) without running out of memory. In the future we might update the code to process larger recordings in smaller chunks.

When you start processing, you should run the processing once while monitoring memory and CPU usage, then adjust the nWorkers/nParallelWorkers parameters to avoid running out of memory during the alignment and summarizing steps.This will greatly improve performance and avoid errors. 
```


```
Tip:
Reprocessing Data.
To reprocess a dataset from scratch, you must delete the trialTable.mat file in your experiment folder, and select overwriteExisting=true in the parameters for the alignment function (either multiROIRegSLAP2 or stripRegBergamo)
```