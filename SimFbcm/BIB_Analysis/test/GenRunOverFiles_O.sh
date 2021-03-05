#! /bin/bash 
#BASEDIR=/afs/cern.ch/work/g/gauzinge/public/BeamGas/cms_conditioned_oxygen

CNT=0
until [ $CNT -gt 111 ]
do
  echo $CNT
  cmsRun BH_generator_cfg_O.py FileNo=$CNT nEvents=1800 >> tmp_O.txt
  rm -rf tmp_O.txt
  ((CNT++))
  
done


