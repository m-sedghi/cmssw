#! /bin/bash 

CNT=33
cmsRun BH_generator_cfg_C.py FileNo=$CNT nEvents=1800 >> tmp.txt
rm -rf tmp.txt
