In order to run CTPPS simulation in CMSSW_8_1_0 follow these steps in terminal:
~~~~
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_9_2_6
cd CMSSW_9_2_6/src
cmsenv
git cms-merge-topic mackoo13:simulation_9_2_6
scram b -j 4
cd Configuration/Test
cmsRun test.py
~~~~

