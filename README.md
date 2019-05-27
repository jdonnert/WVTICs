# WVT ICs

### *Description* 
- Generate glass-like initial conditions for Smoothed Particle Hydrodynamics
- Relaxation of the particle distribution is done using al algorithm based on Weighted Voronoi Tesselations (Diehl et al. 2012)
- Additional particle reshuffling can be enabled to improve over- and undersampled maxima/minima
- A full suite of (analytical) testproblems is included

### *Installation*
Tailor the Makefile (or CMakeLists.txt) to your machine, and run
`make` or `cmake`

### *Usage*
Adjust the parameterfile `ics.par` and run the code. See the `examples` folder for example parameterfiles.
`OMP_NUM_THREADS=4 ./WVTICs ics.par`

### *Contributing*
Contributions are *welcome* and will be *credited*

### *Credits*
Please cite:
- [Donnert 2017, MNRAS](http://adsabs.harvard.edu/abs/2017MNRAS.471.4587D "Toycluster Paper")
- [Arth et al. 2017, MNRAS](http://adsabs.harvard.edu/cgi-bin/basic_connect?qsearch=%5EArth%202017 "WVTICs Paper") 

### *License*
GPL-3.0

### *Examples and post-processing*
- For Python
   - See the Python files in the `scripts` folder

- For IDL
  - Do we want to add example scripts?


(c) Donnert, Arth, Steinwandel, Halbesma 2016-7
