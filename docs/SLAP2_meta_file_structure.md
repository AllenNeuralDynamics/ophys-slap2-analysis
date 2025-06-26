# SLAP2 .meta File Structure Documentation

## Overview

SLAP2 `.meta` files are MATLAB v7.3 HDF5 files that contain acquisition metadata and hardware configuration parameters for SLAP2 experiments. These files are generated during data acquisition and store essential information about the experimental setup, hardware settings, and acquisition parameters.

## File Format

- **Format**: MATLAB v7.3 MAT-file (HDF5 schema)
- **Platform**: PCWIN64
- **Extension**: `.meta`
- **Example**: `acquisition_20250508_145020_DMD1.meta`

## File Structure Hierarchy

```
Root Level
├── #refs# (Group) - Internal MATLAB object references
├── #subsystem# (Group) - MATLAB subsystem information
├── AcquisitionContainer (Group) - Main acquisition configuration
├── machineConfiguration (Group) - Hardware configuration details
└── [Direct Parameters] - Core acquisition and hardware settings
```

## Root Level Parameters

### Core Acquisition Parameters

| Parameter | Type | Description | Example Value |
|-----------|------|-------------|---------------|
| `acqDuration_s` | float64 | Total acquisition duration in seconds | `99999999.` |
| `acqDurationCycles` | float64 (2×1) | Duration in acquisition cycles | `[4.2949673e+09, 4.2949673e+09]` |
| `acquisitionPathIdx` | float64 | Index of the acquisition path | `1.` |
| `acquisitionPathName` | uint16 array | Name of acquisition path (encoded) | `"Path1"` |
| `linePeriod_s` | float64 | Time period per scan line in seconds | `9.35745477e-05` |
| `samplesPerLine` | float64 | Number of samples per scan line | `9261.` |
| `fpgaTimeReference` | uint32 (1×6) | FPGA timing reference array | `[3707764736, 2, 1, 1, 1, 1]` |

### Hardware Settings

| Parameter | Type | Description | Example Value |
|-----------|------|-------------|---------------|
| `aomActive` | uint8 | AOM (Acousto-Optic Modulator) active state | `1` (enabled) |
| `aomVoltage` | float64 | AOM voltage setting | `0.` |
| `channelsSave` | float64 | Number of channels to save | `1.` |
| `enableStack` | float64 | Z-stack acquisition enabled | `0.` (disabled) |
| `remoteFocusPosition_um` | float32 | Remote focus position in micrometers | `56.75` |

### DMD (Digital Micromirror Device) Parameters

| Parameter | Type | Description | Example Value |
|-----------|------|-------------|---------------|
| `dmdPixelsPerColumn` | float64 | Number of DMD pixels per column | `800.` |
| `dmdPixelsPerRow` | float64 | Number of DMD pixels per row | `1280.` |

## AcquisitionContainer Group

Contains detailed acquisition planning and scanning parameters:

### Sub-groups:
- **AcquisitionPlan**: Acquisition sequence planning
- **DmdPatternSequence**: DMD pattern configurations
- **ParsePlan**: Data parsing configuration
  - `acqParsePlan`: ROI and pixel mapping
  - `isSimpleRaster`: Raster scan mode flag
  - `lineRateHz`: Line scan rate (e.g., `10686.66699219`)
  - `linesPerCycle`: Lines per acquisition cycle
  - `pixPerLine`: Pixels per line
  - `rasterOffsetXY`: Raster offset coordinates
  - `rasterSizeXY`: Raster scan dimensions
- **ROIs**: Region of Interest definitions
- **ScannerParameters**: Scanner configuration
  - `fovSize`: Field of view dimensions
  - `lineShear`: Line shear correction values
  - `pixelDilationXY`: Pixel dilation parameters

## machineConfiguration Group

Contains comprehensive hardware configuration with three main arrays:
- **configuration**: Hardware component configurations
- **instanceClass**: Component class definitions  
- **instanceName**: Component instance names

Each array contains 18 object references corresponding to different hardware components including:
- Acquisition paths
- Laser systems
- DMD devices
- Scanner components
- Detection systems
- Motion control hardware

## Hardware Component Categories

### Imaging Components
- **Acquisition Paths**: `Path1`, `Path2` (multi-channel imaging)
- **DMD Systems**: `DMD1`, `DMD2` (spatial light modulation)
- **Remote Focus**: Piezo-based axial positioning

### Laser and Optics
- **AOM Control**: Acousto-optic modulation
- **Laser Clock**: Timing synchronization
- **Pockels Cell**: Fast laser modulation
- **Power Control**: Laser power management

### Scanning Systems
- **Polygonal Scanner**: High-speed line scanning
- **Scanner Control**: Speed and synchronization
- **Motion Detection**: Real-time motion tracking

### Data Acquisition
- **DAQ Systems**: Multi-channel data acquisition
- **Timing Control**: Trigger and clock management
- **Digital I/O**: Control signal routing

## Key Acquisition Parameters Summary

```yaml
Timing:
  - acqDuration_s: Total acquisition time
  - linePeriod_s: Line scan period
  - samplesPerLine: Samples per scan line

Spatial:
  - dmdPixelsPerRow: 1280 pixels
  - dmdPixelsPerColumn: 800 pixels  
  - remoteFocusPosition_um: Z-position

Hardware:
  - aomVoltage: Laser modulation
  - channelsSave: Active channels
  - acquisitionPathName: Detection path
```

## Usage Notes

1. **File Access**: Requires HDF5 readers (h5py, scipy.io with HDF5 support)
2. **Encoding**: String parameters stored as uint16 arrays
3. **References**: Complex objects stored in `#refs#` group with reference pointers
4. **Compatibility**: MATLAB v7.3 format ensures cross-platform compatibility
5. **Validation**: Critical for experiment reproducibility and data provenance

## Related Files

- **Summary Files**: Processed data with experiment results (`*_Summary.mat`)
- **Raw Data**: Acquisition data files (various formats)
- **Session Files**: Experiment session metadata (`session.json`)

This documentation provides the essential structure and parameters needed to understand SLAP2 `.meta` files for experiment analysis and data processing workflows.
