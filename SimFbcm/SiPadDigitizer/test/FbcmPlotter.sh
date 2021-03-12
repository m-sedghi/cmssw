#! /bin/bash 

echo "$(dirname "$0")";
cd "$(dirname "$0")";

export SCRAM_ARCH=slc7_amd64_gcc900;
eval `scramv1 runtime -sh`;

./FbcmPlotter.py $@;
