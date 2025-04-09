# monox_fit

This is the recipe for the Run3 mono-x fit code. For Run2 legacy code, check the tag "legacy_run2"

## Setup combine

Start by setting up [combine v10 in CMSSW 14](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#combine-v10-recommended-version), 
and then also set up [combineHarvester](http://cms-analysis.github.io/CombineHarvester/index.html):

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12

cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v10.0.1
scramv1 b clean; scramv1 b # always make a clean build

cd $CMSSW_BASE/src
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
scram b -j4
```

## Create the workspace

Make sure to source your environment before running anything:
```
source setup_env.sh 
```
To build the workspace and combined model, use the script `build_workspace.py`. This script:

1.	Generates a RooWorkspace from histograms in a ROOT file.
2.	Builds a combined model compatible with CMS Combine.
3.	Creates INFO.txt with input/output checksums and Git metadata.
4.	Symlinks a Makefile into the output folder for producing datacards.


To run it, execute
```bash
python3 build_workspace.py -d <input_dir> -a <analysis> -y <year> -v <variable> [-f <root_folder>] [-t <tag>]
```

The script takes these optional arguments:

| Argument         | Description                                                  | Default         |
|------------------|--------------------------------------------------------------|-----------------|
| `-d`, `--dir`     | Path to directory containing the input ROOT file             | *(required)*    |
| `-a`, `--analysis` | Analysis name: `vbf`, `monojet`, `monov`, etc.              | `"vbf"`         |
| `-y`, `--year`    | Year of the dataset: `2017`, `2018`, or `Run3`               | `"2017"`        |
| `-v`, `--variable`| Observable variable name: `mjj`, `met`, etc.                 | `"mjj"`         |
| `-f`, `--folder`  | Folder inside the ROOT file to look for histograms           | auto-detected   |
| `-t`, `--tag`     | Custom output tag (used in output folder name)               | today’s date    |


Example:

```bash
python3 build_workspace.py \
  -d /path/to/inputs \
  -a vbf \
  -y 2018 \
  -v mjj \
  -t 2025_03_31
```

This will:

- Read `/path/to/inputs/limit_vbf.root`
- Create a workspace and model for category `vbf_2018`
- Output files in:  
  `$FIT_FRAMEWORK_PATH/vbf/2018/2025_03_31/root/`

-

```
vbf/
└── 2018/
    └── 2025_03_31/
        ├── root/
        │   ├── ws_vbf.root                # The RooWorkspace
        │   ├── combined_model_vbf.root    # The Combine-ready model
        │   ├── INFO.txt                   # MD5 + Git info
        │   └── Makefile                   # Symlink for Combine datacards
```