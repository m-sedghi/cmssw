#! /bin/bash 
#BASEDIR=/afs/cern.ch/work/g/gauzinge/public/BeamGas/cms_conditioned_hydrogen

CNT=0
until [ $CNT -gt 100 ]
do
  echo $CNT
  cmsRun BH_generator_cfg_H.py FileNo=$CNT nEvents=1800 >> tmp_H.txt
  rm -rf tmp_H.txt
  ((CNT++))
  
done


