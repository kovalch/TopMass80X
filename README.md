Code for top (mass) analyses at Universität Hamburg
------------------------------------------------------

To get started for run-2 on the NAF:


cmsrel CMSSW_7_6_3_patch2

cd CMSSW_7_6_3_patch2/src

cmsenv

export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/$USER/.cmsgit-cache

git cms-init

git clone -b leptonjets_run2 ssh://git@gitlab.cern.ch:7999/stadie/uhh-top-mass.git TopMass

git clone --single-branch --branch UHH_Run2 https://git.cern.ch/reps/TopAnalysis.git

git remote add origin git@github.com:stadie/cmssw.git

git fetch origin

cp TopMass/Configuration/sparse-checkout .git/info/sparse-checkout

git checkout

scram b -j 12


----------------------------------------------------------
To test it on data and MC:

cmsRun TopMass/Configuration/analysis/analyzeTopHypotheses_cfg.py  mcversion=data

cmsRun TopMass/Configuration/analysis/analyzeTopHypotheses_cfg.py