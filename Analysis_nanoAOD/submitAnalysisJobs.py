import os
import sys
import re
import datetime


mX_min = 400
mX_max = 4000
mX_step = 200
mY_step = 200


failedMassPoints = []
#failedMassPoints = [[1600, 500], [3000, 260], [3000, 300], [3000, 500], [3000, 700], [3000, 900], [3000, 1100], [3000, 1300]]

from optparse import OptionParser

usage = "Usage: python %prog [options] \n\nExample: python %prog -o condor_analysis_jobs -s higgsTripletAnalysis.py [--withNu]"

# input parameters
parser = OptionParser(usage=usage)

parser.add_option('-o', '--output', action='store',
                  dest='output',
                  help='Output folder name to which a time stamp is appended at execution time (This parameter is mandatory)')

parser.add_option('-s', '--script', action='store',
                  dest='script',
                  help='Python script (This parameter is mandatory)')

parser.add_option("--withNu", action="store_true",
                  dest="withNu",
                  help="Include neutrinos in GenJets",
                  default=False)

(options, args) = parser.parse_args()

# make sure all necessary input parameters are provided
if not (options.output and options.script):
    print 'Mandatory parameters missing'
    print ''
    parser.print_help()
    sys.exit(1)


def processMassPoint(mX, mY, folder, script):
    print ('---------------------------------------')
    print ('(mX, mY) = (%i, %i)' % (mX, mY))
    print ('---------------------------------------')

    # create jdl file
    jdl_path = os.path.join(folder, ('TRSM_XToHY_6b_M3_%i_M2_%i.jdl' % (mX, mY)))
    jdl_file = open(jdl_path,'w')
    jdl_content = jdl_template
    jdl_content = re.sub('DUMMY_OUTPUTDIR',folder,jdl_content)
    jdl_content = re.sub('M3MASS',str(mX),jdl_content)
    jdl_content = re.sub('M2MASS',str(mY),jdl_content)
    jdl_content = re.sub('DUMMY_SCRIPT',script,jdl_content)
    jdl_file.write(jdl_content)
    jdl_file.close()

    # submit Condor job
    os.system('condor_submit ' + jdl_path)


jdl_template = """universe = vanilla
Notification = never
Executable = DUMMY_OUTPUTDIR/analysisJob.sh
Output = DUMMY_OUTPUTDIR/TRSM_XToHY_6b_M3_M3MASS_M2_M2MASS.stdout
Error = DUMMY_OUTPUTDIR/TRSM_XToHY_6b_M3_M3MASS_M2_M2MASS.stderr
Log = DUMMY_OUTPUTDIR/TRSM_XToHY_6b_M3_M3MASS_M2_M2MASS.log
Arguments = M3MASS M2MASS DUMMY_OUTPUTDIR DUMMY_SCRIPT
Queue
"""


bash_template = """#!/bin/bash

M3MASS=$1
M2MASS=$2
OUTPUTDIR=$3
SCRIPT=$4

START_TIME=`/bin/date`
echo "Started at $START_TIME"
echo ""

cd ${_CONDOR_SCRATCH_DIR}

source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project CMSSW CMSSW_10_6_0
cd CMSSW_10_6_0/src/
eval `scram runtime -sh`

python ${OUTPUTDIR}/${SCRIPT} --mX=${M3MASS} --mY=${M2MASS} DUMMY_NEUTRINO

mv -v HISTOGRAMS_*.root ${OUTPUTDIR}
exitcode=$?

echo ""
END_TIME=`/bin/date`
echo "Finished at $END_TIME"
exit $exitcode
"""

# make condor job folder
condor_folder = os.path.join(os.getcwd(), options.output + '_' + datetime.datetime.now().strftime("%Y%m%d%H%M%S"))
os.system('mkdir -p ' + condor_folder)

# copy Python script
os.system('cp -p %s %s' % (options.script, condor_folder))

# create Bash script
bash_path = os.path.join(condor_folder,'analysisJob.sh')
bash_script = open(bash_path,'w')
bash_content = bash_template
bash_content = re.sub('DUMMY_NEUTRINO',('--withNu' if options.withNu else ''),bash_content)
bash_script.write(bash_content)
bash_script.close()
os.system('chmod +x ' + bash_path)


if len(failedMassPoints):
    for mp in failedMassPoints:
        processMassPoint(mp[0],mp[1],condor_folder,options.script)
else:
    # loop over mass points
    for mX in range(mX_min, mX_max + mX_step, mX_step):
        for mY in sorted(list(set([260,mX-140])) + range(300, mX-125, mY_step)):
            processMassPoint(mX,mY,condor_folder,options.script)

print ('---------------------------------------')
