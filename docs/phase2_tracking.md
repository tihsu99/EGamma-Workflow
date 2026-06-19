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
```
