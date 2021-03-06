#! /bin/bash 
BASEDIR=/afs/cern.ch/work/m/msedghi/public/BeamInducedBackgrdFbcm/bibSIM/BeamGasOxygen

for FILE in $BASEDIR/* ; do 
echo $FILE; 

cmsRun BIB_DIGI_bh_not_as_2ndSource_SelfMixed_cfg.py filein=$FILE
#break

done

