echo "sourcing LCG stack"
source /cvmfs/sft.cern.ch/lcg/views/LCG_102rc1/x86_64-centos7-gcc11-opt/setup.sh

SFDIR=./data/scaleFactors
if [ ! -d $SFDIR ]; then
    echo "get muon efficiency scale factor repo"
    git clone https://gitlab.cern.ch/cms-muonPOG/muonefficiencies.git ./data/scaleFactors
fi
