#!/bin/bash
# Navigate to working directory
cd /work/jdriesch/CMS_service/rochester/new_rochester/rochester_shire_source

source env.sh
python rochester_rdf.py -4 --process $1
