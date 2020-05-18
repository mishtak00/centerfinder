# centerfinder

Clone this repository by running:
```
git clone https://github.com/mishtak00/centerfinder
```

Create, activate and configure a virtual environment that's ready to host centerfinder by running the following from the terminal while inside the centerfinder directory:
```
python -m venv virtualenv
source virtualenv/bin/activate
pip install -r requirements.txt
```

To run the voting procedure, call the --kernel_radius or -r argument followed by the desired test radius:
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143
```

To apply a vote cut to the centers grid, call the --vote_threshold or -t argument followed by desired vote nr:
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -t 170
```
The default threshold is 0 (no cutting by vote number, every centers grid is preserved).

To subtract the background from the galaxy density grid, call the --background_subtract or -b argument:
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -t 170 -b
```
The default voting procedure won't apply background subtraction unless the user requests with the above.

Add the --params_file or -p argument to change the default file from which the cosmological parameters are loaded to a new file whose name is given as argument:
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -t 170 -b -p params_cmassdr8.json
```

Add the --save or -s argument to save outputs in a new file called "out" (created automatically inside current directory):
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -t 170 -b -p params_cmassdr9.json -s
```

Add the --verbose or -v argument to have sanity checks and feedback on the running process printed to standard output as the program runs in exchange for *slightly* degraded performance (some of these checks are sums over large grids):
```
python cfdriver.py mock_cmassDR9_north_3001.fits -r 143 -t 170 -b -p params_cmassdr9.json -s -v
```
