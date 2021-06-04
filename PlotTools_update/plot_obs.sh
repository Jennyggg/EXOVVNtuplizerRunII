#!/bin/bash
#mass=('verylow' 'low' 'medium' 'high' 'veryhigh' 'inclu')
mass=('veryhigh' 'inclu')
lumi=('inclu')
obs=('sph' 'thrust' 'ndisplace' 'ntrk' 'pt' 'masspertrk')
#obs=('pt' 'masspertrk')
mkdir -p plots
for massregion in "${mass[@]}"
do
  for lumiregion in "${lumi[@]}"
  do
    for Obs in "${obs[@]}"
    do
      root -l -b -q plot_obs_data_MC_lumislice.cpp'("/pnfs/psi.ch/cms/trivcat/store/user/jinw/DY1JetsToLL/CRAB/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/210520_001515/0000/flatTuple_compare.root","/pnfs/psi.ch/cms/trivcat/store/user/jinw/ZeroBias2018/run316187/flatTuple_run316187_v4.root","MCDYWeight_run316187_collision.root","plots/ZeroBias_DYMC_ratio_WColl","DY1JetsToLLPileup","ZeroBias","nPU_DY.root","'$Obs'","'$massregion'","'$lumiregion'",true,true,false)'
    done
  done
done
