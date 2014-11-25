To compile the ntuplizer:

Make a 72X release. Please log into a machine running SL6
$ cmsrel CMSSW_7_2_0
$ cd CMSSW_7_2_0/src
$ cmsenv
checkout the repoistory
$ git clone https://github.com/bchiarito/ntuplizer Analysis
build
$ scram b

To use the ntuplizer:

$ cd Analysis/MiniAnalyzer/python
Here you must run add_to_miniAOD.py to get extra jet collections (CA8).  Make sure to change the input file to your miniAOD source.
$ cmsRun add_to_miniAOD.py
Then you can run the ntuplizer.
$ cmsRun b2g_miniAod_to_EDMtuple.py
