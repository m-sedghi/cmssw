#! /bin/bash 
BASEDIR=/afs/cern.ch/work/m/msedghi/public/bibGeneratorOutput/BeamGasOxygen

for FILE in $BASEDIR/* ; do 
echo $FILE; 

cmsRun BIB_SIM_MS_cfg.py filein=$FILE
#break

done

