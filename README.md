# Quick Setup
## Set up this tools
Please note that the `EGamma-Workflow` do not depends on the `cmssw`, so no need to do `cmsenv` in this sense. By default, we run in `lxplus9`.
```
mkdir EGM-Validation
cd EGM-Validation
git clone https://github.com/tihsu99/EGamma-Workflow.git
```
## Install `TnP` ToolKit.
```
# Note, no need of cmsenv
# To utilizing the jupyter notebook mode, please login with localhost
ssh USERNAME@lxplus.cern.ch -L8787:localhost:8787 # 8787 could be changed

# [First time] Install & Create Clean Virtual Environment for the run.

git clone -b 2024_Studies git@github.com:saumyaphor4252/egamma-tnp.git
cd egamma-tnp/
python3 -m venv egmtnpenv
source egmtnpenv/bin/activate
pip install jupyter
pip install ipython
pip install ipykernel
ipython kernel install --user --name=egmtnpenv
python -m ipykernel install --user --name=egmtnpenv
pip install . --no-cache-dir

# [After installation] Source env
source egmtnpenv/bin/activate
voms-proxy-init --voms cms --valid 100:00

# Initial the jupyternotebook
jupyter lab --no-browser --port [PORT NUMBER]
```
There should be a link like `http://localhost:8xxx/?token=...` displayed in the output at this point, paste that into your browser. You should see a jupyter notebook with a directory listing.

# HLT Setup
## Set up cmssw env
```
cd ../
cmsrel [CMSSW VERSION]
cmsenv
# Optional: some custmized functions are in our `python` directory, could moved under to cmssw dir for running
cp python/XXX.py [cmssw dir]/src/.
# Build your HLTConfiguration etc.
...
```

# Structure
## Final structure
```
EGM-Validation/
├── EGamma-Workflow/
│   ├── submit_condor.py
│   ├── config/
│   │   └── *.yaml
│   ├── python/
│   │   └── (useful functions)
│   └── egamma_tnp/        # installed (official toolkit)
└── CMSSW_X_X_X/           # HLT workspace
```

# Comparison plot
## Yaml Config
The pipelines are full controled by the yaml files. This is the format it should follow. Several examples could be find in `config`.
```yaml
cmssw-dir: [cmssw src path]
eos-dir: [output storage place]

project:
  [VersionTagName-1]:
    dataset:
      - [DAS path]
    GT: [GT]
    process:
      # Those commands will be run one-by-one on the fly. Only the specified "storefile" will be captured from the condor machine to your eos-dir
      # Just example here.
      - cmsDriver.py Phase2 -s L1P2GT,HLT:75e33 --processName=HLTX --conditions auto:phase2_realistic_T33 --geometry ExtendedRun4D110 --era Phase2C17I13M9 --eventcontent FEVTDEBUGHLT --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,HLTrigger/Configuration/customizeHLTforEGamma.customiseEGammaMenuDev  --inputCommands='keep *, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT' --mc --no_exec   --python_filename hlt_phase2_mod.py
      - cmsRun hlt_phase2_mod.py
      - python3 /afs/cern.ch/user/t/tihsu/EGM_HLT_JIRA/CMSSW_16_0_0_pre3/src/makeNtuples_Phase2.py --input-dir "./" -o validation.root -n 10000
      - python3 /afs/cern.ch/user/t/tihsu/EGM_HLT_JIRA/CMSSW_16_0_0_pre3/src/makeNtuples_Phase2_v2.py Phase2_L1P2GT_HLT.root -o detailedvalidation.root
      - cmsDriver.py step2_RECO_MINI -s RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --conditions auto:phase2_realistic_T33 --datatier MINIAODSIM --eventcontent MINIAODSIM --geometry ExtendedRun4D110 --era Phase2C17I13M9 --filein file:Phase2_L1P2GT_HLT.root --fileout file:stepMINI.root --hltProcess HLTX --mc --no_exec --python_filename hlt_ministep.py
      - cmsRun hlt_ministep.py
      - python3 /afs/cern.ch/user/t/tihsu/EGM_HLT_JIRA/CMSSW_16_0_0_pre3/src/lastFilterEffHist_Upgrade_MINIAOD.py stepMINI.root -o tnpNtupler.root

    storefile:
      - validation.root
      - detailedvalidation.root
      - tnpNtupler.root
```
## Produce condor jobs
This create the running script without submission
```
python3 submit_condor.py --n [N proceesed per file] --config config/phase2_tracking.yaml --proxy [userproxy] --farm [Farm]
```

## Submission
Submit the condor and you should have results in the `eos-dir`.
```
condor_submit [Farm_Dir]/condor_jobs.sub
```

## Summary plot [PHASE 2]
Several python plotting functions are stored in `plot_func`.
- Trigger Efficiency: `plotHist.py`.
- Variable Distribution: `plot_Phase2_EGM_Variables_v*.py`

## Produce Trigger Efficiency [run3]
Follow the [instruction](#install-tnp-toolkit) in `TnP` ToolKit, and the jupyter notebook example could be found in `notebook` directory.
