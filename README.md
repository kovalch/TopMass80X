Code for top (mass) analyses at Universität Hamburg
------------------------------------------------------

To get started for run-2 on the NAF:


cmsrel CMSSW_7_4_15

cd CMSSW_7_4_15/src

cmsenv

export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/$USER/.cmsgit-cache

git cms-init

git clone -b leptonjets_run2 ssh://git@gitlab.cern.ch:7999/stadie/uhh-top-mass.git TopMass

git clone --single-branch --branch UHH_Run2 https://git.cern.ch/reps/TopAnalysis.git

git remote add git@github.com:stadie/cmg-cmssw.git

git fetch origin

cp TopMass/Configuration/sparse-checkout .git/info/sparse-checkout

git checkout -b CMGTools-from-CMSSW_7_4_12 origin/CMGTools-from-CMSSW_7_4_12

scram b -j 12


----------------------------------------------------------
To test it on data and MC:

cmsRun TopMass/Configuration/analysis/analyzeTopHypotheses_cfg.py  mcversion=data

cmsRun TopMass/Configuration/analysis/analyzeTopHypotheses_cfg.py