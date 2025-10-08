# StandaloneUnfolding Input Requirements

This document describes the required inputs for using the `StandaloneUnfolding` class in the xsec_analyzer framework.

## Overview

The `StandaloneUnfolding` class provides a simplified interface to perform unfolding without running the full analysis pipeline. It requires four TMatrixD objects and one text file as inputs.

## Required Inputs

### 1. data_signal (TMatrixD)

**Description**: Background-subtracted data event counts in reconstructed (reco) bins

**Dimensions**: Column vector (N_reco × 1)

**Calculation**:
```
data_signal = BNB_data - (EXT + MC_background_prediction)
```

**Prepared by**:
- `SystematicsCalculator::get_measured_events()` → `MeasuredEvents.reco_signal_`

**Location in code**:
- `src/utils/SystematicsCalculator.cxx:1495-1538`
- Formula at lines 1528-1529: `ordinary_data - ext_plus_mc_bkgd`

**Physical meaning**: The measured signal after removing beam-off backgrounds (EXT) and MC-predicted beam-correlated backgrounds

---

### 2. data_covmat (TMatrixD)

**Description**: Total covariance matrix on the background-subtracted data measurement

**Dimensions**: Square matrix (N_reco × N_reco)

**Contents**:
- Statistical uncertainties on data
- Systematic uncertainties:
  - Flux prediction variations
  - Cross-section model variations
  - Detector response variations
- Correlations between all reco bins

**Prepared by**:
- `SystematicsCalculator::get_measured_events()` → `MeasuredEvents.cov_matrix_`

**Location in code**:
- `src/utils/SystematicsCalculator.cxx:1506-1514`
- Extracted from `get_covariances()->at("total")`

**Physical meaning**: Full uncertainty on the measurement including all sources of systematic error and their bin-to-bin correlations

---

### 3. smearcept (TMatrixD)

**Description**: Smearceptance matrix (detector response matrix)

**Dimensions**: Rectangular matrix (N_reco × N_true)

**Formula**:
```
smearcept(r,t) = N(true=t, reco=r) / N(true=t)
```
Where:
- N(true=t, reco=r) = events with true bin t reconstructed in reco bin r
- N(true=t) = total events in true bin t

**Prepared by**:
- `SystematicsCalculator::get_cv_smearceptance_matrix()`

**Location in code**:
- `src/utils/SystematicsCalculator.cxx:1389-1414`
- Uses CV universe's 2D histogram divided by true bin counts

**Physical meaning**:
- Element (r,t) represents the probability that a true event in bin t is reconstructed in bin r
- Encodes both detector efficiency and bin migration effects
- Used to predict reco-space distribution from true-space distribution: `reco = smearcept * true`

---

### 4. prior_true_signal (TMatrixD)

**Description**: Prior expectation for true signal event counts

**Dimensions**: Column vector (N_true × 1)

**Contents**: Central-value MC prediction for signal in each true bin

**Prepared by**:
- `SystematicsCalculator::get_cv_true_signal()`

**Location in code**:
- `src/utils/SystematicsCalculator.cxx:1420-1436`
- Extracted from CV universe's `hist_true_` histogram

**Physical meaning**:
- MC-based prior used by iterative unfolding algorithms
- For D'Agostini: starting point for iterations
- For Wiener-SVD: used in regularization calculation

---

### 5. BlocksFile (Text file)

**Description**: Defines bin grouping for blockwise unfolding

**Purpose**:
When the analysis phase space contains disconnected regions (e.g., separate 1D distributions for muon momentum and muon cos(θ)), bins are grouped into "blocks". Each block is unfolded independently to avoid introducing spurious correlations between physically disconnected observables.

**Format**:
```
<num_true_bins>
<true_bin_index_0> <block_index_0>
<true_bin_index_1> <block_index_1>
<true_bin_index_2> <block_index_2>
...
<num_reco_bins>
<reco_bin_index_0> <block_index_0>
<reco_bin_index_1> <block_index_1>
<reco_bin_index_2> <block_index_2>
...
```

**Example**:
```
60
0 0
1 0
2 0
...
29 0
30 1
31 1
...
59 1
57
0 0
1 0
...
28 0
29 1
30 1
...
56 1
```
This example has 60 true bins and 57 reco bins split into 2 blocks.

**Location in code**:
- Read in `src/utils/Unfolder.cxx:274-304` (blockwise_unfold method)
- Parsing logic at lines 287-299

**Physical meaning**:
- Block index groups bins measuring the same type of observable
- Block 0 might be muon momentum bins, Block 1 might be muon angle bins
- Prevents the unfolding from incorrectly correlating momentum and angle measurements

---

## Where These Matrices Are Created in the Full Pipeline

### Standard Analysis Flow

1. **ProcessNTuples** (`bin/ProcessNTuples`):
   - Reads PeLEE ntuples
   - Applies selection cuts
   - Creates `stv_tree` ROOT files with analysis variables

2. **univmake** (`bin/univmake`):
   - Reads `stv_tree` files for all samples
   - Creates systematic universe histograms
   - Output: `Universes.root` containing:
     - 2D histograms (true vs reco) for migration matrices
     - 1D histograms (true, reco separately)
     - Multiple universes for each systematic variation

3. **SystematicsCalculator** class:
   - Constructor: `src/utils/SystematicsCalculator.cxx:54-193`
   - Loads universe histograms from `Universes.root`
   - Computes covariance matrices via `get_covariances()` method
   - Provides access to the four required matrices

4. **CrossSectionExtractor** class:
   - Wrapper used in `src/app/Unfolder.C:34`
   - Internally creates `SystematicsCalculator`
   - Method `get_unfolded_events()` (line 366) calls unfolding

### Key Methods for Matrix Extraction

| Matrix | Method | File Location |
|--------|--------|---------------|
| data_signal | `SystematicsCalculator::get_measured_events()` | `src/utils/SystematicsCalculator.cxx:1495` |
| data_covmat | `SystematicsCalculator::get_measured_events()` | `src/utils/SystematicsCalculator.cxx:1495` |
| smearcept | `SystematicsCalculator::get_cv_smearceptance_matrix()` | `src/utils/SystematicsCalculator.cxx:1389` |
| prior_true_signal | `SystematicsCalculator::get_cv_true_signal()` | `src/utils/SystematicsCalculator.cxx:1420` |

---

## How to Prepare Matrices for Standalone Unfolding

**Important**: The standard analysis pipeline does **not** automatically save these matrices in the format required by `StandaloneUnfolding`. You must create them manually.

### Step 1: Create a Matrix Extraction Script

Create a custom C++ script (e.g., `extract_unfolding_inputs.C`) that loads your universes file and saves the required matrices:

```cpp
// extract_unfolding_inputs.C
#include "XSecAnalyzer/SystematicsCalculator.hh"
#include "XSecAnalyzer/MCC9SystematicsCalculator.hh"
#include "TFile.h"
#include "TMatrixD.h"

void extract_unfolding_inputs(
    const std::string& univ_file,
    const std::string& syst_config,
    const std::string& output_file)
{
    // Create SystematicsCalculator with your universe file
    MCC9SystematicsCalculator syst(univ_file, syst_config);

    // Extract the four required matrices
    auto meas = syst.get_measured_events();
    auto smearcept = syst.get_cv_smearceptance_matrix();
    auto prior = syst.get_cv_true_signal();

    // Save to ROOT file
    TFile out(output_file.c_str(), "recreate");
    out.WriteObject(meas.reco_signal_.get(), "data_signal");
    out.WriteObject(meas.cov_matrix_.get(), "data_covmat");
    out.WriteObject(smearcept.get(), "smearcept");
    out.WriteObject(prior.get(), "prior_true_signal");
    out.Close();

    std::cout << "Saved unfolding inputs to " << output_file << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage: extract_unfolding_inputs UNIV_FILE SYST_CONFIG OUTPUT_FILE\n";
        return 1;
    }

    extract_unfolding_inputs(argv[1], argv[2], argv[3]);
    return 0;
}
```

### Step 2: Create the BlocksFile

The blocks information comes from your binning scheme. You can extract it programmatically:

```cpp
// extract_blocks.C
#include "XSecAnalyzer/SystematicsCalculator.hh"
#include <fstream>

void extract_blocks(
    const std::string& univ_file,
    const std::string& syst_config,
    const std::string& blocks_file)
{
    MCC9SystematicsCalculator syst(univ_file, syst_config);

    const auto& true_bins = syst.true_bins_;
    const auto& reco_bins = syst.reco_bins_;

    std::ofstream out(blocks_file);

    // Write true bins
    out << true_bins.size() << "\n";
    for (size_t i = 0; i < true_bins.size(); ++i) {
        if (true_bins[i].type_ == kSignalTrueBin) {
            out << i << " " << true_bins[i].block_index_ << "\n";
        }
    }

    // Write reco bins
    size_t num_ordinary_reco = 0;
    for (const auto& rb : reco_bins) {
        if (rb.type_ == kOrdinaryRecoBin) ++num_ordinary_reco;
    }

    out << num_ordinary_reco << "\n";
    for (size_t i = 0; i < reco_bins.size(); ++i) {
        if (reco_bins[i].type_ == kOrdinaryRecoBin) {
            out << i << " " << reco_bins[i].block_index_ << "\n";
        }
    }

    out.close();
    std::cout << "Saved blocks to " << blocks_file << std::endl;
}
```

### Step 3: Create StandaloneUnfolding Configuration

Create a config file (e.g., `standalone_unfold.conf`):

```
InputFile unfolding_inputs.root
OutputFile unfolded_result.root
BlocksFile blocks.txt
Unfold WienerSVD 1 second-deriv
```

Or for D'Agostini unfolding:
```
InputFile unfolding_inputs.root
OutputFile unfolded_result.root
BlocksFile blocks.txt
Unfold DAgostini iter 4
```

### Step 4: Run Standalone Unfolding

```bash
# Extract the matrices
./extract_unfolding_inputs Universes.root configs/systcalc.conf unfolding_inputs.root

# Extract the blocks
./extract_blocks Universes.root configs/systcalc.conf blocks.txt

# Run the unfolding
bin/StandaloneUnfold standalone_unfold.conf
```

---

## Configuration File Format

The StandaloneUnfolding configuration file supports the following commands:

### InputFile
```
InputFile <path_to_root_file>
```
ROOT file containing the four required TMatrixD objects.

### OutputFile
```
OutputFile <path_to_output_file>
```
ROOT file where unfolding results will be saved.

### BlocksFile
```
BlocksFile <path_to_blocks_file>
```
Text file defining bin blocks (format described above).

### Unfold

**Wiener-SVD unfolding:**
```
Unfold WienerSVD <use_filter> <reg_type>
```
- `<use_filter>`: 0 or 1 (whether to use Wiener filter)
- `<reg_type>`: identity, first-deriv, or second-deriv

**D'Agostini iterative unfolding:**
```
Unfold DAgostini iter <num_iterations>
```
Or:
```
Unfold DAgostini fm <figure_of_merit>
```
- `iter`: specify fixed number of iterations
- `fm`: iterate until figure of merit reaches threshold

---

## Output

The `StandaloneUnfolding::run_unfolding()` method saves the following objects to the output ROOT file:

| Object Name | Type | Description |
|-------------|------|-------------|
| `unfolded_signal` | TMatrixD | Unfolded event counts in true bins (N_true × 1) |
| `cov_matrix` | TMatrixD | Covariance matrix on unfolded signal (N_true × N_true) |
| `unfolding_matrix` | TMatrixD | Unfolding matrix M such that true = M * reco (N_true × N_reco) |
| `add_smear_matrix` | TMatrixD | Additional smearing matrix from regularization (N_true × N_true) |

See `src/utils/StandaloneUnfolding.cxx:174-203` for implementation details.

---

## References

### Key Source Files

- `include/XSecAnalyzer/StandaloneUnfolding.hh` - Header file
- `src/utils/StandaloneUnfolding.cxx` - Implementation
- `src/utils/Unfolder.cxx` - Base unfolding algorithms
- `src/utils/SystematicsCalculator.cxx` - Matrix preparation
- `src/app/standalone_unfold.C` - Main executable

### Related Documentation

- See `CLAUDE.md` for overall framework architecture
- See `configs/tutorial_bin_config.txt` for binning scheme examples
- See `configs/xsec_config.txt` for CrossSectionExtractor configuration examples
