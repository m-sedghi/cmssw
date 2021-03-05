#! /bin/bash 
#BASEDIR=/afs/cern.ch/work/g/gauzinge/public/BeamGas/cms_conditioned_carbon

CNT=0
until [ $CNT -gt 111 ]
do
  echo $CNT
  cmsRun BH_generator_cfg_C.py FileNo=$CNT nEvents=1800 >> tmp_C.txt
  rm -rf tmp_C.txt
  ((CNT++))
  
done


