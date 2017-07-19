In order to run CTPPS simulation in CMSSW_8_1_0 follow these steps in terminal:
~~~~
source /afs/cern.ch/cms/cmsset_default.sh
cmsrel CMSSW_8_1_0_pre8
cd CMSSW_8_1_0_pre8/src
cmsenv
git cms-merge-topic ziaber:simulation
scram b -j 4
cd Configuration/Test
cmsRun test.py
~~~~

