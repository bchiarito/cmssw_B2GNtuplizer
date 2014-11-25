
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

Open b2g_miniAod_to_EDMtuple and change the input file and output name to what you want. Then run the ntuplizer.

$ cmsRun b2g_miniAod_to_EDMtuple.py
