# Quick Setup
## Set up this tools
```
mkdir EGM-Validation
cd EGM-Validation
git clone https://github.com/tihsu99/EGamma-Workflow.git
```
## Set up cmssw env
```
cd ../
cmsrel [CMSSW VERSION]
cmsenv
```

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

## Produce Trigger Efficiency
Will put the instruction here soon. Quick example could be found in `notebook` 
