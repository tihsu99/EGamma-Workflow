import os
import yaml
import subprocess
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Prepare Condor jobs for HLT rerun")
parser.add_argument("--jobFlavour", type=str, default="tomorrow",
                    help="Condor jobflavor (e.g., espresso, longlunch, workday, tomorrow)")
parser.add_argument("--n", type=int, default=10, help="Number of DAS files per shell script")
parser.add_argument("--farm", type=str, default="Farm",
                    help="farm directory")
parser.add_argument("--nEvent", type=int, default=-1)
parser.add_argument("--config", type=str)
parser.add_argument("--proxy", type=str, default=None)
args = parser.parse_args()
jobflavor = args.jobFlavour

# Load YAML
with open(args.config) as f:
    cfg = yaml.safe_load(f)

cmssw_dir = cfg['cmssw-dir']
eos_dir = cfg['eos-dir']

shell_scripts = []

def get_das_files(dataset):

  das_command = f'dasgoclient -query="file dataset={dataset}"'
  files = os.popen(das_command).read().strip().split("\n")
  accessible = []

  for file in files:
    file_path = f"/eos/cms/{file}"

    if os.path.exists(file_path):
        accessible.append(f"file:{file_path}")
    
  return accessible



for project_name, project_cfg in cfg['project'].items():
    project_eos_dir = os.path.join(eos_dir, project_name)
    os.makedirs(project_eos_dir, exist_ok=True)

    farm_dir = os.path.join(args.farm, project_name)
    os.makedirs(farm_dir, exist_ok=True)
    
    all_input_files = []
    for ds in project_cfg['dataset']:
        all_input_files.extend(get_das_files(ds))
    
    start_idx = 0
    end_idx = min(args.n, len(all_input_files))
    bunch_idx = 0
    while (end_idx < (len(all_input_files))):
        script_name = os.path.join(farm_dir, f"run_{project_name}_{bunch_idx}.sh")
        workspace = f"$TMPDIR/Job_{project_name}_{bunch_idx}"
        print(f"creating {script_name}")
        with open(script_name, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("export X509_USER_PROXY=$1\n")
    #        f.write("voms-proxy-info -all\n")
    #        f.write("voms-proxy-info -all -file $1\n")
            f.write(f"WORKDIR={workspace}\n")
            f.write(f"mkdir -p $WORKDIR\n")
            f.write("echo 'Using TMPDIR=' $WORKDIR\n")
            f.write(f"cd {cmssw_dir}\n")
            f.write("eval `scramv1 runtime -sh`\n\n")
            f.write("cd $WORKDIR\n")
            for file_idx in range(start_idx, end_idx):
                for i, cmd in enumerate(project_cfg['process']):
                    cmd_tmp = cmd.replace("file:", f"file:$WORKDIR/")
                    if "TnPTreeProducer" in cmd_tmp:
                        cmd_tmp = cmd_tmp.replace('outputFile=file:', 'outputFile=')
                    if i == 0:
                       cmd_tmp += f" --filein {all_input_files[file_idx]} "
                    if 'cmsDriver' in cmd_tmp:
                       cmd_tmp += f" --customise_commands 'process.maxEvents.input=cms.untracked.int32({args.nEvent})' "
                    f.write(f"echo 'Running step {i+1}'\n")
                    f.write(f"{cmd_tmp}\n")
        
                for storefile in project_cfg.get('storefile', []):
                    eos_path = os.path.join(project_eos_dir, storefile.replace(".root", f"_{file_idx}.root"))
                    f.write(f"xrdcp $WORKDIR/{storefile} root://eosuser.cern.ch/{eos_path}\n")
                f.write("rm $WORKDIR/*.root\n")
    
        os.chmod(script_name, 0o755)
        shell_scripts.append(script_name)
        end_idx += args.n
        start_idx += args.n
        bunch_idx += 1


sub_file = f"{args.farm}/condor_jobs.sub"
condor_str = "executable = $(filename)\n"
if args.proxy is not None:
    condor_str += f"Proxy_path = {args.proxy}\n"
    condor_str += "arguments = $(Proxy_path)\n"
condor_str += "output = $Fp(filename)$(filename)hlt.stdout\n"
condor_str += "error = $Fp(filename)$(filename)hlt.stderr\n"
condor_str += "log = $Fp(filename)$(filename)hlt.log\n"
condor_str += f'+JobFlavour = "{args.jobFlavour}"\n'
condor_str += f"queue filename matching ({args.farm}/*/*.sh)"
condor_file = open(sub_file, "w")
condor_file.write(condor_str)

