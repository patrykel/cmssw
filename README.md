In order to run CTPPS simulation in CMSSW_8_1_0 follow these steps in terminal:
~~~~
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_8_1_0
cd CMSSW_8_1_0/src
cmsenv
git cms-merge-topic mackoo13:simulation_8_1_0
scram b -j 8
cd Configuration/Test
cmsRun test.py
~~~~
The code is based on official CMS 8_1_0 release and Marcin ziaber's fork (ziaber/simulation)
