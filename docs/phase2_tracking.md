# Phase2 Tracking Setup
This is the setup for the Phase2 tracking performance study, which is used to perform the release validation.
## Setup
```bash
# under the working area of your choice
cd WORKDIR
cmsrel CMSSW_16_0_0_pre3
cd CMSSW_16_0_0_pre3/src
cmsenv
git cms-addpkg HLTrigger/Configuration
cp [EGamma_WORKFLOW_DIR]/python/customizeHLTforEGamma.py [WORKDIR]/HLTrigger/Configuration/python/.
scram b -j 8
```
## Prepare the workflow config
```aiignore
# Replace the path in config/phase2_tracking.yaml to your own path
# Some python files are used to produce the ntuples, which we have provided in the python directory i.e. makeNtuples_Phase2.py etc.
python3 submit_condor.py --config config/phase2_tracking.yaml
# Then you can submit the condor jobs to the farm or run it locally
```
