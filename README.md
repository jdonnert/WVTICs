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

### *Example MpsFractions*
| Problem           |   1E4    |   1E5    |  1E6     |
| :---------------: | :------: | :------: | :------: |
| Constant Density  | 1.7      | 0.8      | 0.35     |
|    TopHat         | 3.5      | 1.7      | 0.8      |
|    SawTooth       | 0.15     | 0.08     | 0.04     |
|    Sinus          | 1.3      |  0.6     | 0.3      |
|    Mach 2         | 0.4      | 0.19     | x        |
|    Mach 3         | 0.4      | 0.19     | x        |
|    Mach 4         | 0.4      | 0.19     | x        |
|    Sod Shock      | 5.0e-6   | 0.19     | x        |
|    Sedov Blast    | 1.7      | 0.8      | 0.35     |
| Kelvin-Helmholtz  | 15.0     | 9.0      | 12.0     |
|  Keplarian Ring   | 14.0     | 6.0      | x        |
|    Cold Blob      | 0.5      | 0.23     | 0.1      |
| Evrard Collapse   | 1.7      | 0.8      | 0.35     |
| Zel’dovic Pancake | 1.7      | 0.8      | 0.35     |
|      Box          | 14.0     | 6.5      | 3.2      |
|   Gresho-Vortex   | 11.5     | 5.5      | x        |
| Exponential Disk  | 80.0     | x        | x        |
|  Boss-Bodenheimer | 0.6      | 0.3      | 0.15     |
|   Beta Profile    | 1.7      | 0.8      | 0.35     |
|    Fast Rotor     | 12.0     | 6.0      | 3.5      |
|    Strong Blast   | 12.0     | 6.0      | 3.5      |
| Orszang-Tang Vort.| 12.0     | 6.0      | 3.5      |
| Lin. Alfvén Wave  | 0.2      | x        | x        |


### *Contributing*
Contributions are *welcome* and will be *credited*

### *Credits*
Please cite:
- [Donnert 2014, MNRAS](http://adsabs.harvard.edu/abs/2014MNRAS.438.1971D "Toycluster Paper")
- [Arth et al. 2017, MNRAS](http://adsabs.harvard.edu/cgi-bin/basic_connect?qsearch=%5EArth%202017 "WVTICs Paper")

### *License*
GPL-3.0

### *Examples and post-processing*
- For Python
   - See the Python files in the `scripts` folder

- For IDL
  - Do we want to add example scripts?


(c) Donnert, Arth, Steinwandel, Halbesma 2016-7
