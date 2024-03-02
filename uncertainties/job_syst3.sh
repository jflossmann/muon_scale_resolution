#!/bin/bash
# Navigate to working directory
cd /work/jdriesch/CMS_service/rochester/new_rochester/rochester_shire_source

source env.sh
python rochester_rdf_syst.py -1 -2 -3 -4 --process $1 --syst 3
