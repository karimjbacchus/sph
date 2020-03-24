# ACSE-4-SPH

[Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) (SPH) is a meshless
method for solving the Navier-Stokes equation, in which fluid properties are stored on Lagrangian fluid particles (i.e. on
particles which move with the fluid flow). The particles interact to generate values across the entire fluid domain through
continuous smoothing kernels. 

As the SPH method is meshless and Lagrangian, it is ideal for solving problems involving fluid flow with interfaces and free 
surfaces. This tool implements the SPH method in C++ to solve wave generation in a lock-release/dam-break problem.

### Compilation/Installation Guide:

Clone the repository from Github by either

* Using the command line:

```
git clone https://github.com/acse-2019/acse-4-sph-doon.git
```
* Downloading the repository as a .zip

Next, within the command prompt enter
```
pip install -r requirements.txt
```
To compile using command line arguments for a g++ compiler from within the directory, run:
```
g++ src/*.cpp -o a
```
where ```*``` is the file name. Then run:
```
./a
```
to execute the program.

### User instructions:

#### Running the simulation
Before running the simulation, the user can set the values of the initial particle spacing ```dx```, the speed of sound ```c0```, the initial density ```rho0```, the factor for the smoothing length ```h_fac```, and the time step ```dt```. These values can be set in ```SPH_2D.cpp```in the function ```SPH_main::set_values```. The user can also set the duration of the simulation (in seconds) by setting the value of ```maxt``` in ```SPH_Snippet.cpp``` or ```SPH_Snippet_Euler.cpp```.

The user has two time-stepping methods available: a Forward Euler method implemented in ```SPH_Snippet_Euler.cpp``` and a Predictor-Corrector method implemented in ```SPH_Snippet.cpp```.

The user can also choose between different initial configurations, as detailed in the documentation, by choosing between ```domain.place_points```,```domain.place_points_conf1```, ```domain.place_points_conf2```,```domain.place_points_conf3```, ```domain.place_points_conf4```, ```domain.place_points_conf5``` in ```SPH_Snippet.cpp```.

Once the cpp file has been compiled and executed, the output files are stored in the folder ```output```, with each file representing an iteration. Note that the file ```output_-1.vtp``` is used for testing, and it is not an output of the simulation.

#### Viewing the simulation
The results and its analysis can be viewed in the Jupyter Notebook ```ACSE-4.3-htmlview.ipynb```. 

This notebook contains five Python functions which allow the user to read the output files, get the physical quantities of the simulation -namely the velocity and pressure-, create and fill grids for each physical quantity, and measure the properties at the crest of the wave. Each of these functions are described in detail in the Sphinx documentation as well as in ``UserGuide_ACSE-4.3-htmlview.pdf``. To run all the functions at once, use ```ensemble(no, tol)```. 

Once the user has input their chosen conditions, they can then plot the wave itself at a specified timestep and animate it over a period of time.

For a full **convergence analysis** of the simulator, open the Jupyter notebook ```Convergence_analysis.ipynb```.


### Documentation:

The code includes [Sphinx](https://www.sphinx-doc.org) documentation. On systems with Sphinx installed, this can be built by running

```
python -m sphinx docs docs/html
```

then viewing the `index.html` file in the `docs/html` directory in your browser.


### Testing


The tool includes tests, which you can use to check its operation on your system. With the code compiled, these can be run 
with

```
python run_tests.py
```

Alternatively, on the command line enter:
```
make runtests
```
to run all the tests.


### Credits

Team Doon:
Karim Bacchus, Luca Cilio, Hanran Ji, Nadya Mere, Tianzong Yu.
