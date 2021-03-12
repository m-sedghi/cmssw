eval `scramv1 runtime -sh`;
source /cvmfs/cms.cern.ch/common/crab-setup.sh;
for pu in 200 140 100 50 10 0p5 1 1p5
do 
    crab submit -c crabNTuple_cfg.py General.requestName=FBCMnTuplePU$pu
done
