# CenterFinder

### Contents
* [Setup](#setup)
* [I/O](#io)
* [Run](#run)
* [Exit](#exit)

### Setup
CenterFinder is written purely in Python 3. Please make sure to have it installed in your system before proceeding with the rest of these instructions. This readme assumes that the `python` keyword has been aliased to `python3`. This readme contains code that's meant for Unix-based terminals (Mac, Linux...). Windows analogues may be *slightly* different. 

Clone this repository by running the following in the terminal:
```
git clone https://github.com/mishtak00/centerfinder
```

*Optional, __recommended__:* Create and activate a virtual environment to host CenterFinder by running the following from the terminal while inside the centerfinder directory:
```
python -m venv cf_env
source cf_env/bin/activate
```

You can choose any package manager for python to install CenterFinder's requirements in the virtual environment (or directly to your system, if you didn't follow the *optional recommendation* above). I have chosen pip here as an example:
```
pip install -r requirements.txt
```

### I/O
CenterFinder expects its input data in Astropy's FITS format (python recarrays). It reads in a 3 or 4-column [TableHDU](https://docs.astropy.org/en/stable/io/fits/api/tables.html#) depending on whether the input has a 4th column whose values are to serve as weights. The principal output of CenterFinder is a 4-column .fits table that contains location and weight information on each found baryon acoustic oscillation center. Each row is a 4-tuple of the right ascension, declination, redshift and weight of a center.

Additionally, the program expects a JSON file named 'params.json'. The 'params_cmassdr9.json' file has been given here as a template parameters file. This file is expected to contain select cosmological parameters as well as a crucial class variable for CenterFinder: grid_spacing. Grid spacing is the side length of each cubic cell in the big grids that represent histograms in real space in CenterFinder. The latter is suggested to be kept between 5 and 10 h-1Mpc during testing stages, selecting values under 5 only when more refined analysis is needed. CenterFinder's runtime is inversely proportional to the cube of this variable, i.e. linear in grid volume.

CenterFinder's main routine accepts weighted or unweighted input and deals with both cases accordingly. Run CenterFinder on weighted input through the -w or --weighted_input argument. In this case, the program expects a 4-column table. Otherwise, feed CenterFinder a 3-column input, omit the -w argument and it will treat each input data point read as having unit weight. The user can also feed a weighted catalog (4-column) to CenterFinder and omit -w to ignore the weights from the catalog and treat each data point as unit-weighted. A 3-column BinTableHDU containing CMASS DR9 mock data has been provided here as a test case.


### Run
To run the voting procedure on unweighted input, just call the --kernel_radius or -r argument followed by the desired test radius. This refers to the *"2.2.1 Galaxy Count"* portion of the paper.
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143
```
To run the voting procedure on weighted input, call the -w argument with the -r argument followed by the desired test radius. This refers to the *"2.2.2 Weighted Count"* portion of the paper.
```
python cfdriver.py mock_cmassDR9_north_3001.fits -w -r 143
```

To subtract the expected grid from the galaxy density grid, call the --density_contrast or -c argument like below. This refers to the *"2.2.3 Density Contrast"* portion of the paper. The default voting procedure won't apply background subtraction unless the user requests it.
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -c
```

To subtract the average weighted density from the galaxy density grid, call the --overdensity or -o argument like below. This refers to the *"2.2.4 Over-Density"* portion of the paper. The default voting procedure won't apply the overdensity correction unless the user requests it. **Note:** The -c and -o arguments are mutually exclusive. Only density contrast or over-density can be used at a time to generate the density field, not both.
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -o
```

To apply a vote cut to the centers grid, call the --vote_threshold or -t argument followed by desired vote number. The default threshold is 0 (no cutting by vote number, every bin in the centers grid is preserved). **Note:** An original, uncut centers grid will be saved to the outputs folder if the user specifies the -s argument (see below).
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -c -t 170
```

Add the --params_file or -p argument to change the default file from which the cosmological parameters are loaded to a new file whose name is given as argument:
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -c -t 170 -p params_cmassdr8.json
```

Add the --save or -s argument to save outputs in the "out_INPUT_FILENAME" folder (created automatically inside current directory):
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -c -t 170 -p params_cmassdr9.json -s
```

Add the --verbose or -v argument to have sanity checks and feedback on the running process printed to standard output as the program runs:
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -c -t 170 -p params_cmassdr9.json -s -v
```

### Exit
To deactivate the virtual environment once the job has finished:
```
deactivate
```
