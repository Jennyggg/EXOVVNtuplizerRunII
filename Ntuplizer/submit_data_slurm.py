import os

TEMPLATE='''#!/bin/bash
# 
#SBATCH -p standard
#SBATCH --account=t3
#SBATCH --job-name=data_process
#SBATCH --mem=3G                     # memory 3GB (per job)
#SBATCH --time 12:00:00


echo HOME: $HOME 
echo USER: $USER 
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME

# each worker node has local /scratch space to be used during job run
mkdir -p /scratch/$USER/${{SLURM_JOB_ID}}
export TMPDIR=/scratch/$USER/${{SLURM_JOB_ID}}

#########################################################
# Please fill in your code here
cmsenv
sed -i "s|TMPDIR|$TMPDIR|g" {CONFIG}
cmsRun  {CONFIG}  RunPeriod="Run2018C"
xrdcp $TMPDIR/*root {TARGET}
rm {CONFIG}
#########################################################
# cleaning of temporary working dir after job is completed:
rm  -rf /scratch/$USER/${{SLURM_JOB_ID}}'''

target_loc = "root://t3se01.psi.ch:1094//store/user/jinw/flatTuple_v6/"

data_loc = "/pnfs/psi.ch/cms/trivcat/store/data/Run2018C/ZeroBias/MINIAOD/12Nov2019_UL2018_LowPU-v1/"
for subdir in os.listdir(data_loc):
  for fname in os.listdir(data_loc+subdir):
    with open("config_generic_opt_skimmed_crab.py",'r') as f_read:
      filedata = f_read.read()
    filedata = filedata.replace("OUTPUTFILE","TMPDIR/flatTuple_ZeroBias_2018lowPU_hltv6_"+subdir+"_"+fname)
    filedata = filedata.replace("FOLDER",subdir)
    filedata = filedata.replace("FILENAME",fname)
    fname_write="config_"+subdir+"_"+fname.replace(".root",".py")
    with open(fname_write,'w') as f_write:
      f_write.write(filedata)
    script_name = "batch_"+subdir+"_"+fname.replace(".root",".sh")
    with open(script_name,'w') as script_run:
      script_run.write(TEMPLATE.format(
                        CONFIG = fname_write,
                        TARGET = target_loc,
                      ))
    os.system("sbatch "+script_name)

