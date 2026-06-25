# postProcess_v2

A C++17 post-processing tool that converts raw singles data from **PenRed** Monte Carlo PET simulations into coincidence events and list-mode (LM) output files, ready for image reconstruction.

## Overview

PenRed produces per-detector singles files. This tool reads those files, applies energy and time blurring, runs a sliding-window coincidence sorter, projects hits onto the detector face, and writes a structured binary list-mode stream together with optional histograms (normalization maps, randoms maps, etc.).

```
PenRed singles (binary) ──▶  penred2Coincidence  ──▶  LM file (.lm)
                                                   ──▶  Histograms (normalization, randoms…)
```

## Features

- **Multiple coincidence methods** (selectable via config):
  | ID | Method |
  |----|--------|
  | 0 | Highest-energy single wins (`TAKE_WINNER_IF_ALL_ARE_GOODS`) |
  | 1 | Closest in time (`TAKE_CLOSEST`) |
  | 2 | Same simulation history (`TAKE_SAME_HISTORY`) |
  | 3 | Same history, restricted to 511 keV singles, and closest in energy (`TAKE_SAME_HISTORY_511`) |
  | 4 | Multiplexing — uses prompt gamma to resolve triples (`MULTIPLEXING`) |

- **Multiple projection methods** to map 3-D hit positions onto the detector pixel grid:
  | ID | Method |
  |----|--------|
  | 0 | None (use raw coordinates) |
  | 1 | Fit in detector face |
  | 2 | Extend detector face |

- **Gaussian blurring** of energy, position, and time to simulate realistic detector resolution.
- **Multiple handling** — configurable discard of triples / all multiples.
- **Randoms identification** via history tracking.
- **LM header** (`LMHeader`) carrying acquisition metadata (isotope, activity, timing, geometry…).
- Optional generation of **normalization histograms** (LOR-index maps).

## Requirements

| Tool | Minimum version |
|------|-----------------|
| C++ compiler (GCC / Clang) | C++17 support |
| CMake | 3.16 |

No external libraries are required.

## Building

```bash
# Create and enter the build directory
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

For a debug build with AddressSanitizer and UBSan:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j$(nproc)
```

The resulting binary is `build/penred2Coincidence`.

## Usage

```bash
./build/penred2Coincidence <configFile> <input_prefix>
```

| Argument | Description |
|----------|-------------|
| `configFile` | Path to the configuration file (see format below) |
| `input_prefix` | Path prefix for the PenRed singles files (one file per detector module) |

The tool will look for files named `<input_prefix>0`, `<input_prefix>1`, … (one per detector module).

### Example

```bash
./build/penred2Coincidence config/template_config /data/sim/singles_
```

## Configuration File

Configuration is a plain key-value text file (comments start with `#`). A fully annotated template is provided at [`config/template_config`](config/template_config).

### Key Parameters

```ini
# --- General Simulation ---
Counts          = 100000000000   # Total number of simulated primaries
IsotopeName     = FDG
IsotopeHalfLife = 6586.26        # seconds
AcquisitionTime = 1              # seconds

# --- Geometry ---
DetectorSizeX           = 4.80000E+00  # cm
DetectorSizeY           = 4.80000E+00  # cm
DistanceBetweenDetectors= 1.17000E+01  # cm
ModulesPerRing          = 8
NumberOfRings           = 3
DistanceBetweenRings    = 3.00000E-01  # cm
CrystalDepth            = 1.0          # cm

# --- Resolution & Windows ---
EnergyResolution  = 6.07e-6   # fractional FWHM at 511 keV
PositionResolution= 0.00002   # mm
TimeResolution    = 500.0e-12 # s (FWHM)
EnergyWindow      = 50.0      # % around 511 keV
CoincidenceWindow = 5.0e-9    # s

# --- Bins & Pixels ---
ModuleReconstructionBinsX = 300
ModuleReconstructionBinsY = 300
DetectorPixelsX           = 32
DetectorPixelsY           = 32

# --- Post-Processing Options ---
OutputFormat      = 1   # 0→LM all events, 1→LM coincidences only
CoincidenceMethod = 2   # see table above
ProjectionMethod  = 0   # see table above
PairListFilename  = ../config/detectorList_oldInsert.txt
UseLogicalDetectors = 0
DiscardMultiples    = 0 # 0→keep, 1→discard triples, 2→discard all multiples

# --- Output Flags ---
GenerateHistogram = 0   # 0→No, 1→LOR-index histogram
SaveWeight        = 0   # match your PenRed tally settings
SaveMetadata      = 0
```

### Pair List Files

The `PairListFilename` points to a text file that lists valid detector-module pairs (opposite-facing modules that form lines of response). Pre-built lists are provided in `config/`:

| File | Description |
|------|-------------|
| `detectorList_oldInsert.txt` | Old insert geometry (small) |
| `detectorListInsert12x6.txt` | 12×6 insert geometry |
| `pairListInsertHw.txt` | Hardware insert pair list |

## Project Structure

```
postProcess_v2/
├── CMakeLists.txt
├── config/
│   ├── template_config          # Annotated configuration template
│   ├── detectorListInsert12x6.txt
│   ├── detectorList_oldInsert.txt
│   └── pairListInsertHw.txt
└── src/
    ├── main.cpp                 # Entry point & main processing loop
    ├── common/
    │   ├── types.hh             # Core structs: SingleEvent, CoincidenceEvent, LMHeader, DetectorPair
    │   ├── constants.hh         # Physical and numerical constants
    │   └── math_utils.hh        # Matrix operations, coordinate transforms
    ├── config/
    │   ├── config.hh            # SimConfig struct
    │   └── config.cpp           # Config file parser
    ├── coincidence/
    │   ├── coincidence_engine.hh
    │   └── coincidence_engine.cpp  # Sliding-window coincidence sorter
    ├── geometry/
    │   ├── geometry.hh/.cpp     # Ring/module geometry builder
    │   └── projection.hh/.cpp   # Hit-to-detector-pixel projection
    ├── histogram/
    │   ├── histograms.hh/.cpp   # Accumulation of normalization/randoms maps
    │   └── histogram_writer.hh/.cpp
    ├── io/
    │   ├── singles_reader.hh/.cpp  # Reads & merges per-module PenRed singles files
    │   └── pair_list.hh/.cpp       # Loads detector pair lists
    └── output/
        ├── lm_writer.hh
        └── lm_writer.cpp        # Binary LM file writer
```

## Output

- **List-mode file** (binary): sequence of `LMHeader` + `CoincidenceEvent` records.  
  Each `CoincidenceEvent` contains: timestamp, pair ID, pixel coordinates (x, y) for both detectors, energies, amount, and a gate flag (1 = triple/prompt-gamma event).
- **Histogram files** (if `GenerateHistogram ≠ 0`): written as `modules_*.bin` files in the working directory.

## Statistics Printed at End

```
Number of coincidences: …
Number of singles     : …
Number of randoms     : …
Number of triples     : …
Number of quadruples  : …
Number of 5+ multiples: …
```
