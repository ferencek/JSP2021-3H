# To print out decay kinematics
python XtoYHto3H_decay.py --mX 3000 --mY 3000


# To print the mass points and get the total number
python printMassPoints.py


# To plot mass points
python plotMassPoints.py


# To inspect the generated event(s)
cmsRun inspectEvents_cfg.py maxEvents=1 skipEvents=0 inputFiles=file:/STORE/ferencek/TRSM_XToHY_6b/2017/13TeV/GEN/TRSM_XToHY_6b_M3_3000_M2_300_GEN.root


# To run analysis
python higgsTripletAnalysis.py --maxEvents=-1 --reportEvery=100 --mX=3000 --mY=300 [--withNu]

The last argument --withNu is optional and will process GenJets with included neutrinos.
The mass point can also be defined as, e.g. --massPoint=M3_3000_M2_300 or --massPoint=M3_2500_M2_300_BPf


# To run analysis for all mass points on Condor
python submitAnalysisJobs.py -o condor_analysis_jobs -s higgsTripletAnalysis.py [--withNu]


# To produce analysis plots
python examplePoints_1D.py
python examplePoints_2D.py
python makePlots.py


# To check for missing files and get a list of failed mass points
python checkAnalysisOutput.py path/to/output/files/
