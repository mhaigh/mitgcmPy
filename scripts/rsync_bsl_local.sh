#!/bin/bash
# A variation of script on Archer2. Use to this get data from BAS linux to local machine.

BSL=michai@bslcenb.nerc-bas.ac.uk
BSLDIR=${BSL}:/data/oceans_output/shelf/michai/mitgcm
HOMEDIR=/home/michai/Documents/data #/Users/mh115/Documents/BAS/data
DATAPATH=/PISOMIP_001/run/

NC=1
if [ $NC -eq 1 ] ; then
  # rsync the data/meta files.
  rsync -azvl ${BSLDIR}${DATAPATH}/stateTheta* ${HOMEDIR}${DATAPATH}
fi

rsync -avzL ${BSLDIR}${DATAPATH}/*.meta ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/*.data ${HOMEDIR}${DATAPATH}

# Without SSH keys set up, this will ask for password for each rsync command.

exit

# rsync peripheral files
rsync -avzL ${BSLDIR}${DATAPATH}/hFac* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/Depth* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/DRF* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/DRC* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/DXG* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/DYG* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/RAC* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/RC* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/RF* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/XC* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/YC* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/XG* ${HOMEDIR}${DATAPATH}
rsync -avzL ${BSLDIR}${DATAPATH}/YG* ${HOMEDIR}${DATAPATH}
#rsync -avzL ${BSLDIR}${DATAPATH}/stdout_* ${HOMEDIR}${DATAPATH}

NC=1
if [ $NC -eq 1 ] ; then
  # rsync the data/meta files.
  rsync -azvl ${BSLDIR}${DATAPATH}/stateTheta* ${HOMEDIR}${DATAPATH}
fi                                         
