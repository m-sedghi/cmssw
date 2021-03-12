eval `scramv1 runtime -sh`;
source /cvmfs/cms.cern.ch/common/crab-setup.sh;

for aging in 1000 3000 4000
do
    for pu in 200 140 100 50 10 0p5 1 1p5
    do 
	#echo General.requestName=FBCMDIGIAging${aging}PU${pu}
	crab submit -c crabDIGI_cfg.py General.requestName=FBCMDIGIAging${aging}PU${pu}
    done
done
