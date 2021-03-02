#! /bin/bash 
#BASEDIR=/afs/cern.ch/work/g/gauzinge/public/BeamHalo

CNT=0
until [ $CNT -gt 306 ]
do
  echo $CNT
  cmsRun BH_generator_cfg.py FileNo=$CNT nEvents=1600 >> tmp_BH.txt
  rm -rf tmp_BH.txt
  ((CNT++))
  
done


