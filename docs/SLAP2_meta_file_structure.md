# SLAP2 Meta File Structure

This document describes the structure and contents of SLAP2 `.meta` files. These are MATLAB v7.3 HDF5 files that contain acquisition metadata and hardware configuration settings for SLAP2 experiments.

## File Format

- **Format**: MATLAB v7.3 MAT-file (HDF5)
- **Extension**: `.meta`
- **Typical size**: ~800 KB
- **Platform**: PCWIN64

## Overview

The `.meta` file stores acquisition parameters, hardware settings, and configuration information for each SLAP2 experiment. The file is structured as an HDF5 file with nested groups and datasets containing various acquisition and hardware parameters.

## Parameter Reference

The `.meta` file contains the following root-level parameters:

| Parameter Name | Data Type | Description | Example Value |
|-----------|------|-------------|---------------|
| **Acquisition Parameters** | | | |
| `acqDuration_s` | float64 | Total acquisition duration in seconds | `99999999.` |
| `acqDurationCycles` | float64 (2×1) | Duration in acquisition cycles | `[4.2949673e+09, 4.2949673e+09]` |
| `acquisitionPathIdx` | float64 | Index of the acquisition path | `1.` |
| `acquisitionPathName` | uint16 array | Name of acquisition path (encoded) | `"Path1"` |
| `linePeriod_s` | float64 | Time period per scan line in seconds | `9.35745477e-05` |
| `samplesPerLine` | float64 | Number of samples per scan line | `9261.` |
| `fpgaTimeReference` | uint32 (1×6) | FPGA timing reference array | `[3707764736, 2, 1, 1, 1, 1]` |
| **Hardware Control** | | | |
| `aomActive` | uint8 | AOM (Acousto-Optic Modulator) active state | `1` (enabled) |
| `aomVoltage` | float64 | AOM voltage setting | `0.` |
| `channelsSave` | float64 | Number of channels to save | `1.` |
| `enableStack` | float64 | Z-stack acquisition enabled | `0.` (disabled) |
| `remoteFocusPosition_um` | float32 | Remote focus position in micrometers | `56.75` |
| **DMD Parameters** | | | |
| `dmdPixelsPerColumn` | float64 | Number of DMD pixels per column | `800.` |
| `dmdPixelsPerRow` | float64 | Number of DMD pixels per row | `1280.` |
| **Complex Data Structures** | | | |
| `AcquisitionContainer` | Group | Main acquisition configuration and ROI definitions | Complex nested structure |
| `machineConfiguration` | Group | Hardware configuration details | Complex nested structure |

## Notes

- The file contains many additional internal MATLAB references in the `#refs#` and `#subsystem#` groups
- Detailed hardware configuration parameters are stored within the `machineConfiguration` group
- ROI definitions and acquisition plans are contained within the `AcquisitionContainer` group
- All timing and acquisition parameters are precisely stored for experiment reproducibility

## AcquisitionContainer Sub-fields

The `AcquisitionContainer` group contains detailed acquisition configuration:

| Sub-field | Type | Description |
|-----------|------|-------------|
| **Acquisition Planning** | | |
| `AcquisitionPlan` | Dataset (2,) uint64 | Acquisition sequence planning references |
| `DmdPatternSequence` | Dataset (2,) uint64 | DMD pattern sequence references |
| **Parse Configuration** | | |
| `ParsePlan` | Group | Data parsing and ROI configuration |
| `ParsePlan/acqParsePlan` | Group | Acquisition parsing plan with ROI mappings |
| `ParsePlan/isSimpleRaster` | uint8 | Simple raster scan mode flag |
| `ParsePlan/lineRateHz` | float64 | Line scan rate in Hz (e.g., `10686.66699219`) |
| `ParsePlan/linesPerCycle` | uint64 | Lines per acquisition cycle |
| `ParsePlan/linesPerFrame` | uint32 | Lines per frame |
| `ParsePlan/pixPerLine` | uint32 | Pixels per line |
| `ParsePlan/rasterOffsetXY` | uint32 (2×1) | Raster offset coordinates |
| `ParsePlan/rasterSizeXY` | uint32 (2×1) | Raster scan dimensions |
| `ParsePlan/zs` | float32 | Z-position |
| **ROI and Scanner Configuration** | | |
| `ROIs` | Group | Region of interest definitions |
| `ROIs/rois` | object (4×1) | ROI object references |
| `ScannerParameters` | Group | Scanner configuration parameters |
| `ScannerParameters/fovSize` | float64 (2×1) | Field of view size |
| `ScannerParameters/lineRateHz` | float32 | Line scan rate |
| `ScannerParameters/lineShear` | float32 (1×800) | Line shear correction values |
| `ScannerParameters/maxDmdPatterns` | float64 | Maximum DMD patterns |
| `ScannerParameters/pixelDilationXY` | float64 (2×1) | Pixel dilation parameters |

## machineConfiguration Sub-fields

The `machineConfiguration` group contains hardware component configurations:

| Sub-field | Type | Description |
|-----------|------|-------------|
| `configuration` | object (1×18) | Hardware component configuration array |
| `instanceClass` | object (1×18) | Component class definitions |
| `instanceName` | object (1×18) | Component instance names |

### What "18 total" means

Each of the three arrays (`configuration`, `instanceClass`, `instanceName`) contains **18 object references**. These references point to detailed configuration data stored elsewhere in the file (in the `#refs#` group). 

The 18 hardware components configured in this SLAP2 system include:
- **2 Acquisition Paths** (Path1, Path2)
- **2 DMD devices** (DMD1, DMD2) 
- **2 Remote focus piezo axes** (Axis 1, Axis 2)
- **Various laser control components** (AOM, laser clock, Pockels cells)
- **Scanner system** (polygonal scanner with control channels)
- **Detection systems** (photodiodes, digitizers)
- **Motion control hardware** (motor stages)
- **Data acquisition hardware** (DAQ channels)

## Understanding the #refs# Group

The `#refs#` group is a special MATLAB internal structure used in HDF5-based MAT files (v7.3). It serves as a **reference storage system** for complex MATLAB objects.

**For most users and analysis pipelines, you do NOT need to access or document the contents of `#refs#`.**

## Usage

These files can be read using:
- MATLAB: `load('filename.meta', '-mat')`
- Python: `h5py.File('filename.meta', 'r')`
- HDF5 tools: Any HDF5-compatible reader

For analysis pipelines, focus on the root-level parameters listed above, as they contain the essential acquisition metadata needed for data processing.
