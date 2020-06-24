# CenterFinder

Clone this repository by running the following in the terminal:
```
git clone https://github.com/mishtak00/centerfinder
```

Create, activate and configure a virtual environment that's ready to host centerfinder by running the following from the terminal while inside the centerfinder directory:
```
python -m venv virtualenv
source virtualenv/bin/activate
pip install -r requirements.txt
```

CenterFinder expects its input data in Astropy's FITS format (python recarrays). It reads in a 3 or 4-column TableHDU depending on whether the input has a 4th column whose values are to serve as weights.

Additionally, the program expects a JSON file named 'params.json'. The 'params_cmassdr9.json' file has been given here as a template parameters file. This file is expected to contain select cosmological parameters as well as a crucial class variable for CenterFinder: grid_spacing. Grid spacing is the side length of each cubic bin in the big grids that represent histograms in real space in CenterFinder. The latter is suggested to be kept between 5 and 10 h-1Mpc during testing stages, selecting values under 5 only when more refined analysis is needed. CenterFinder's runtime is inversely proportional to the cube of this variable, i.e. linear in grid volume.

CenterFinder's main routine accepts weighted or unweighted input and deals with both cases accordingly. Run CenterFinder on weighted input through the -w or --weighted_input argument. In this case, the program expects a 4-column table. Otherwise, just omit it and CenterFinder will treat each input data point read as having unit weight. A 3-column BinTableHDU containing CMASS DR9 mock data has been provided here as a test case.

To run the voting procedure on unweighted input, just call the --kernel_radius or -r argument followed by the desired test radius. This refers to the "2.2.1 Galaxy Count" portion of the paper.
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143
```
To run the voting procedure on weighted input, call the -w argument with the -r argument followed by the desired test radius. This refers to the "2.2.2 Weighted Count" portion of the paper.
```
python cfdriver.py mock_cmassDR9_north_3001.fits -w -r 143
```

To subtract the expected grid from the galaxy density grid, call the --background_subtract or -b argument like below. This refers to the "2.2.3 Density Contrast" portion of the paper. The default voting procedure won't apply background subtraction unless the user requests it.
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -b
```

To apply a vote cut to the centers grid, call the --vote_threshold or -t argument followed by desired vote nr as shown below. The default threshold is 0 (no cutting by vote number, every centers grid is preserved).
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -b -t 170
```

Add the --params_file or -p argument to change the default file from which the cosmological parameters are loaded to a new file whose name is given as argument:
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -b -t 170 -p params_cmassdr8.json
```

Add the --save or -s argument to save outputs in a new file called "out" (created automatically inside current directory):
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -b -t 170 -p params_cmassdr9.json -s
```

Add the --verbose or -v argument to have sanity checks and feedback on the running process printed to standard output as the program runs:
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -b -t 170 -p params_cmassdr9.json -s -v
```

To deactivate the virtual environment once the job has finished:
```
deactivate
```
