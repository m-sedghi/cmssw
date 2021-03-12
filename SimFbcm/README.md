# Setup workspace for FBCM simulation in CMSSW_11_2_X
Create your CMS_11_2_X and after setting scram environment variables, please go to release base directory.
Copy the setup script to the CMSSW_11_2_X base directory and run the script.

* Please follow something like the following commands:

cd YourWorkDir

cmsrel CMSSW_11_2_0_pre10

cd CMSSW_11_2_0_pre10/src

cmsenv

cd ..

./Setup_Fbcm_CMSSW_11_2_X.sh


* Keep in mind that the Setup_Fbcm_CMSSW_11_2_X script should be run in the CMSSW_11_2_0_pre10 directory
