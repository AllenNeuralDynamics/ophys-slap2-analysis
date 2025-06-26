# SLAP2 Experimental Summary File Structure Documentation

## Overview

SLAP2 Summary files are MATLAB v7.3 HDF5 files containing processed glutamate imaging data and acquisition metadata. These files store the complete experimental results from the SLAP2 glutamate imaging pipeline.

## File Format Details

- **Format**: MATLAB v7.3 (HDF5-based)
- **Typical Size**: 8-10 GB
- **Extension**: `.mat`
- **Naming**: `Summary-YYMMDD-HHMMSS.mat`

## File Structure Hierarchy

SLAP2 Summary files contain multiple top-level structures:

```
Summary-YYMMDD-HHMMSS.mat
├── exptSummary/                    # Main processed experimental data
│   ├── E (cell array)             # Trial data
│   ├── params (struct)            # Processing parameters  
│   ├── trialTable (struct)        # Trial organization
│   └── [other analysis results]
└── #refs#/                         # MATLAB object references (internal)
└── #subsystem#/                    # MATLAB metadata (internal)
```

**Note**: Summary files contain only processed experimental data. Original acquisition parameters (like `acqDuration_s`, `aomVoltage`, etc.) are stored in separate `.meta` files from the original SLAP2 acquisition.

## ExperimentSummary Structure (exptSummary)

The main data structure containing all processed experimental results.

### Top-Level Fields

| Field | Type | Description |
|-------|------|-------------|
| **Core Data Fields** | | |
| `E` | Cell Array (2×155) | Main experimental data - traces and analysis for each FOV×trial |
| `meanIM` | Cell Array (2×1) | Average images across all trials for each FOV |
| `actIM` | Cell Array (2×1) | Activity images showing identified synaptic sites for each FOV |
| `params` | Structure | Processing parameters and settings used |
| **Trial Organization** | | |
| `trialTable` | Structure | Trial metadata, filenames, timing, and organization |
| `Z` | Array (2×1) | Z-positions for each FOV in micrometers |
| **Additional Analysis Data** | | |
| `userROIs` | Cell Array (2×1) | User-defined regions of interest for each FOV |
| `aData` | Cell Array (2×155) | Additional analysis data per trial |
| `peaks` | Cell Array (2×1) | Peak detection results for each FOV |
| **Per-Trial Images** | | |
| `perTrialActIMs` | Cell Array (2×1) | Activity images for each individual trial |
| `perTrialActIMsAligned` | Cell Array (2×1) | Aligned activity images for each trial |
| `perTrialMeanIMs` | Cell Array (2×1) | Mean images for each individual trial |
| `perTrialMeanIMsAligned` | Cell Array (2×1) | Aligned mean images for each trial |
| `perTrialAlignmentOffsets` | Cell Array (2×1) | Alignment offset data for each trial |
| **Metadata** | | |
| `dr` | String Array (94×1) | Directory path information |

## E{} Cell Array Structure (Main Data)

Each `exptSummary.E{fov, trial}` contains processed data for one field of view in one trial.

### Data Fields in Each E{} Cell

| Field | Type | Description |
|-------|------|-------------|
| `dF` | Structure | Baseline-subtracted fluorescence signals |
| `F0` | Array (e.g., 4002×120) | Fluorescence baseline estimates |
| `dFF` | Structure | Fractional fluorescence changes (dF/F0) |
| `ROIs` | Structure | User-defined region of interest signals |
| `footprints` | Array | Spatial footprints of identified synaptic sources |
| `global` | Structure | Global experiment parameters |
| `noiseEst` | Array | Noise estimation data |
| `discardFrames` | Array | Frames marked for exclusion from analysis |

### Signal Processing Subfields

The `dF`, `dFF`, and `ROIs` structures contain subfields with different temporal processing methods:

| Subfield | Description |
|----------|-------------|
| `ls` | Least Squares - no temporal prior or constraints |
| `matchFilt` | Matched Filter - optimal for event detection |
| `nonneg` | Non-negative - mild nonnegativity constraint |
| `spikes` | Spike Detection - sharp event onsets, removes indicator kinetics |
| `denoised` | Denoised - strong temporal prior, smooth signals |

## Parameters Structure (params)

Processing parameters used to generate the data:

| Parameter | Example Value | Description |
|-----------|---------------|-------------|
| `analyzeHz` | 200.0 | Analysis sampling rate (Hz) |
| `alignHz` | 80.0 | Alignment sampling rate (Hz) |
| `sigma_px` | 1.33 | Expected synaptic source size (pixels) |
| `tau_s` | 0.05 | Signal rise time constant (seconds) |
| `tau_full` | 10.0 | Full signal decay time (seconds) |
| `dXY` | 4.0 | Spatial downsampling factor |
| `nParallelWorkers` | 8.0 | Number of parallel processing workers |
| `sz` | [323, 832] | Image dimensions (pixels) |
| `maxSynapseDensity` | 0.01 | Maximum allowed synapse density |
| `baselineWindow_Glu_s` | 4.0 | Baseline estimation window (seconds) |
| `denoiseWindow_s` | 0.25 | Denoising window size (seconds) |
| `motionThresh` | 2.5 | Motion detection threshold |
| `nanThresh` | 0.25 | NaN tolerance threshold |
| `sparseFac` | 0.05 | Sparsity factor for matrix factorization |

## Trial Table Structure (trialTable)

Organization and metadata for all trials:

| Field | Type | Description |
|-------|------|-------------|
| `epoch` | Array (155×1) | Epoch number for each trial |
| `filename` | Cell Array (155×2) | Source filenames for each trial and FOV |
| `firstLine` | Array (155×2) | First scan line number for each trial |
| `lastLine` | Array (155×2) | Last scan line number for each trial |
| `fnAdata` | Cell Array (155×2) | Analysis data filenames |
| `fnRaw` | Cell Array (155×2) | Raw data filenames |
| `fnRegDS` | Cell Array (155×2) | Registered downsampled filenames |
| `trialStartTimeInferred` | Array (155×1) | Inferred trial start times |
| `trialEndTimeFromPC` | Array (155×1) | Trial end times from PC clock |
| `trueTrialIx` | Array (155×1) | True trial indices |
| `refStack` | Cell Array (2×1) | Reference stack information for each FOV |
| `alignParams` | Structure | Alignment parameters used |

### Alignment Parameters (trialTable.alignParams)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `alignHz` | 80.0 | Alignment sampling rate |
| `alpha` | 0.005 | Alignment regularization parameter |
| `clipShift` | 5.0 | Maximum shift clipping value |
| `includeIntegrationROIs` | false | Include integration ROIs in alignment |
| `isReVolt` | false | Whether using ReVolt system |
| `maxshift` | 50.0 | Maximum allowed shift (pixels) |
| `nWorkers` | 16.0 | Number of alignment workers |
| `overwriteExisting` | false | Whether to overwrite existing files |
| `refStackTemplate` | false | Use reference stack as template |

## Additional Processing Parameters (params)

The `params` structure contains processing parameters used to generate the Summary file. Based on actual file analysis:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `activityChannel` | 1.0 | Channel used for activity detection |
| `alignHz` | 80.0 | Alignment sampling rate (Hz) |
| `analyzeHz` | 200.0 | Analysis sampling rate (Hz) |
| `baselineWindow_Ca_s` | 4.0 | Calcium baseline estimation window (seconds) |
| `baselineWindow_Glu_s` | 4.0 | Glutamate baseline estimation window (seconds) |
| `dXY` | 4.0 | Spatial downsampling factor |
| `denoiseWindow_s` | 0.25 | Denoising window size (seconds) |
| `discardInitial_s` | 0.1 | Initial time to discard (seconds) |
| `drawUserRois` | true | Whether user ROIs were drawn |
| `maxSynapseDensity` | 0.01 | Maximum allowed synapse density |
| `microscope` | [array] | Microscope configuration parameters |
| `motionThresh` | 2.5 | Motion detection threshold |
| `nParallelWorkers` | 8.0 | Number of parallel processing workers |
| `nanThresh` | 0.25 | NaN tolerance threshold |
| `nmfIter` | 5.0 | Number of NMF iterations |
| `numChannels` | 1.0 | Number of acquisition channels |
| `sigma_px` | 1.33 | Expected synaptic source size (pixels) |
| `sparseFac` | 0.05 | Sparsity factor for matrix factorization |
| `sz` | [323, 832] | Image dimensions (pixels) |
| `tau_full` | 10.0 | Full signal decay time (seconds) |
| `tau_s` | 0.05 | Signal rise time constant (seconds) |
