Higgs triplet analysis HOW-TO
#################

1. python submitAnalysisJobs.py -o condor_nanoAOD_[mjet/msoftdrop] -s higgsTripletAnalysis_nanoAOD.py [--msoftdrop]

    -submits jobs to condor to run analysis on all mass points and creates a folder with output histograms
    -run once for jet mass, once for soft drop mass (using respective options)

2. wait

3. python makePlots.py --omjet condor_nanoAOD_mjet_TIMESTAMP --omsoftdrop condor_nanoAOD_msoftdrop_TIMESTAMP 

    -uses analysis output from both jet mass and soft drop mass to create plots 


#################
Everything
#################

# To print out decay kinematics
python XtoYHto3H_decay.py --mX 3000 --mY 3000


# To print the mass points and get the total number
python printMassPoints.py


# To plot mass points
python plotMassPoints.py


# To inspect the generated event(s)
cmsRun inspectEvents_cfg.py maxEvents=1 skipEvents=0 inputFiles=file:/STORE/ferencek/TRSM_XToHY_6b/2017/13TeV/NANOAOD/TRSM_XToHY_6b_M3_4000_M2_300_NANOAOD.root


# To run analysis 
python higgsTripletAnalysis_nanoAOD.py --maxEvents=-1 --reportEvery=100 --mX=4000 --mY=300 [--withNu] [--msoftdrop]

# The argument --withNu is optional and will process FatJets with included neutrinos.
# The argument --msoftdrop is optional and will search for higgs candidates with condition on soft drop mass rather than on pure jet mass


# To run analysis for all mass points on Condor where -o is path to output folder; for analysis with softdrop mass use option --msoftdrop; default is pure jet mass 
python submitAnalysisJobs.py -o condor_nanoAOD_m[jet/softdrop] -s higgsTripletAnalysis_nanoAOD.py [--msoftdrop]


# Takes histograms with mjet AND msoftdrop analysis and produces analysis plots into output folder figs; use --omjet and --omsoftdrop to specify the folder which contains analysis output
python makePlots.py --omjet condor_nanoAOD_mjet_TIMESTAMP --omsoftdrop condor_nanoAOD_msoftdrop_TIMESTAMP


# To check for missing files and get a list of failed mass points
python checkAnalysisOutput.py path/to/output/files/

# To produce analysis plots
python examplePoints_1D.py
python examplePoints_2D.py