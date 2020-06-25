#!/bin/bash
#SBATCH --partition=standard -c 1 --mem-per-cpu=60gb -t 2:00:00 --output=out.log
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -t 170 -c -p params_cmassdr9.json -s -v
