# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a cross-section analysis framework for MicroBooNE neutrino physics measurements. It processes PeLEE ntuples from the "searchingfornues" framework to extract differential cross-sections using systematic uncertainty propagation and unfolding techniques. The framework implements selection-based analysis with full support for detector systematics, flux variations, and multi-universe reweighting for precision neutrino interaction measurements.

## Essential Commands

### Environment Setup
```bash
# Source ROOT and environment (select appropriate method for your system)
source setup_xsec_analyzer.sh  # Auto-detects OS and sets up ROOT

# For AlmaLinux 9
source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
spack load gcc@12.2.0 arch=linux-almalinux9-x86_64_v3
spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3

# For Scientific Linux 7
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode v08_00_00_84 -q e17:prof
```

### Build Commands
```bash
# Standard build (optimized)
make

# Debug build (with symbols and -O0)
make debug

# Clean build
make clean && make

# Build creates:
# - lib/libXSecAnalyzer.so (shared library)
# - bin/ProcessNTuples (ntuple processor)
# - bin/univmake (universe histogram maker)
# - bin/Unfolder (unfolding application)
# - bin/SlicePlots (slice histogram plotter)
# - bin/xsroot (ROOT REPL with framework loaded)
# - bin/xsnotebook (Jupyter notebook launcher)
# - bin/AddFakeWeights (NuMI dirt file preprocessor)
# - bin/AddBeamlineGeometryWeights (NuMI beamline systematic weights)
# - bin/UnfolderNuMI (NuMI-specific unfolding)
```

### Running Analysis

```bash
# Interactive ROOT with framework loaded
xsroot

# Process ntuples to create analysis trees
ProcessNTuples input.root file_type selection_name output.root

# Create universe histograms for systematics
univmake config/file_properties.txt config/bin_config.txt output_universes.root

# Produce slice plots with systematics
SlicePlots config/file_properties.txt config/systcalc.conf config/slice_config.txt universes.root output/

# Run unfolding
Unfolder config/xsec_config.txt config/slice_config.txt output/ measurement.root

# Full analysis pipeline
scripts/FullAnalysis.sh
```

### NuMI-Specific Workflows

To enable NuMI mode, set `useNuMI = true` in `include/XSecAnalyzer/Constants.hh` and rebuild.

```bash
# Preprocess NuMI overlay files with beamline geometry weights
AddBeamlineGeometryWeights input.root FHC output_with_weights.root

# Add fake weight maps to dirt files (required for systematics framework)
AddFakeWeights dirt_input.root reference_overlay.root output_dirt.root

# Run NuMI unfolding
UnfolderNuMI config/xsec_config_numi.txt config/slice_config.txt output/ measurement.root
```

## High-Level Architecture

### Core Analysis Pipeline

1. **ProcessNTuples**: First-stage processing
   - Reads PeLEE ntuples from searchingfornues framework
   - Applies selection cuts via `SelectionBase` implementations
   - Computes true and reconstructed observables
   - Outputs simplified ROOT TTree (`stv_tree`) with analysis variables

2. **univmake**: Systematic universe creation
   - Processes `stv_tree` files for all data/MC samples
   - Creates multi-universe histograms for reweightable systematics (flux, cross-section)
   - Handles detector systematic variations (detVar samples)
   - Produces migration matrices (true vs reco binning)

3. **SlicePlots**: Measurement binning
   - Projects multi-dimensional distributions into analysis bins ("slices")
   - Calculates covariance matrices across systematic universes
   - Handles background subtraction and efficiency corrections

4. **Unfolder**: Cross-section extraction
   - Implements Wiener SVD and D'Agostini iterative unfolding
   - Corrects for detector effects using migration matrices
   - Propagates systematic uncertainties through unfolding
   - Outputs differential cross-sections with full covariance

### Selection Framework

All physics selections inherit from `SelectionBase` and implement:
- `define_signal()`: Truth-level signal definition
- `selection()`: Reconstruction-level event selection
- `categorize_event()`: Event classification for backgrounds
- `compute_reco_observables()`: Observable calculation from reco
- `compute_true_observables()`: Observable calculation from truth

Selections are instantiated via `SelectionFactory` and configured through text files.

Examples:
- `CC1mu1p0pi`: CC 1-muon 1-proton 0-pion selection
- `CC1muNp0pi`: CC 1-muon N-proton 0-pion selection
- `NuMICC1e`: NuMI electron neutrino CC selection

### Binning System

**BinScheme classes** (in `src/binning/`) define:
- True-level bin definitions (for MC signal)
- Reco-level bin definitions (for data/selected events)
- Sideband regions (control samples)
- Organized into "blocks" representing different phase space regions

**SliceBinning** projects N-dimensional bins into 1D "slices" for plotting and unfolding.

### Systematic Uncertainty Handling

**Reweightable systematics** (MCC9SystematicsCalculator):
- GENIE cross-section model parameters
- Flux prediction uncertainties
- Handled via multi-universe weights stored in ntuples

**Detector systematics** (separate MC samples):
- Wire response variations
- Space charge effect variations
- Electron lifetime variations
- Handled by comparing detVar samples to detVar CV

### Configuration Files

- `configs/files_to_process.txt`: Input ntuple file list
- `configs/file_properties.txt`: Sample metadata (POT, normalization, type)
- `configs/*_bin_config.txt`: True/reco bin definitions
- `configs/*_slice_config.txt`: Slice variable definitions and binning
- `configs/systcalc.conf`: Systematic uncertainty configuration
- `configs/xsec_config.txt`: Cross-section extraction settings

### Key Data Structures

- `AnalysisEvent`: Container for ntuple branch data
- `Universe`: Single systematic variation with histograms
- `UniverseMaker`: Creates and manages systematic universes
- `SliceHistogram`: Multi-dimensional histogram projected to 1D
- `FiducialVolume`: Detector fiducial volume definition (6 boundaries)

### Technology Stack

- **ROOT**: C++ physics analysis framework (TTree, TH1, TH2, etc.)
- **C++17**: Primary implementation language
- **GNU Make**: Build system
- **Jupyter**: Interactive analysis via `xsnotebook`

### Important Constants

Located in `include/XSecAnalyzer/Constants.hh`:
- Fiducial volume boundaries (FV_X_MIN, FV_X_MAX, etc.)
- Particle PDG codes
- Physics cuts (momentum thresholds, PID scores)
- Particle masses and binding energy
- `useNuMI`: Toggle BNB/NuMI beam mode (requires rebuild)

### Workflow Example

```bash
# 1. Process ntuples for all samples
ProcessNTuples pelee_overlay.root overlay CC1mu1p0pi stv_overlay.root
ProcessNTuples pelee_data.root data CC1mu1p0pi stv_data.root

# 2. Create systematic universes
univmake configs/file_properties.txt configs/cc1mu_bin_config.txt universes.root

# 3. Make slice plots with uncertainties
SlicePlots configs/file_properties.txt configs/systcalc.conf configs/slice_config.txt universes.root plots/

# 4. Extract cross-section
Unfolder configs/xsec_config.txt configs/slice_config.txt plots/ xsec_result.root
```

### Special Considerations

- All executables require `XSEC_ANALYZER_DIR` environment variable (set by `setup_xsec_analyzer.sh`)
- Selections are compiled into `libXSecAnalyzer.so` and registered via `SelectionFactory`
- Systematic weights must follow specific naming conventions (`weight_*_0`, `weight_*_1`, etc.)
- Detector variations use unweighted universes with different detector response
- NuMI mode changes flux calculation, systematic handling, and available selections
