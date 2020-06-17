#!/bin/bash
#SBATCH --partition=standard -c 1 --mem-per-cpu=32gb -t 2:00:00 --output=out.log
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -t 170 -b -p params_cmassdr9.json -s -v
