[![Code Climate](https://codeclimate.com/github/Stefan-Endres/DWPM-Mixture-Model/badges/gpa.svg)](https://codeclimate.com/github/Stefan-Endres/DWPM-Mixture-Model)

# DWPM-Mixture-Model
Phase seperation calculation using the DWPM mixture rule.

## Basic description of files.
1. `main.py` main script to call, the specifications in `config.cfg` is used to 
determine the scripts to be run. See Config files for detail.
2. `nComp.py` contains the core functions used to simulate multicomponent 
equilibrium problems and optimise parameter models. 
3. `pure.py` is used to fit single component Van der Waals EOS temperature 
dependent activity coefficient parameters to the Adachi-Lu and Soave models. 
These parameters are used in the multicomponent functions. The single component 
data is contained in `\Data\pure` and every component simulated in 
multicomponent system requires vapour pressure data or parameters stores in the 
.csv files in `\Data\pure`. 
4. `Van_der_Waals.py` contains the functions the volume root, saturation pressure and state dependent functions which are called in both `pure.py` and `nComp.py` (note that for example, the function `a_T(s, p)` is used in calculating the temperature dependent activity coefficient parameters for both the pure and the mixed phases states (for example to calculated for component 1: `s.c[1] = a_T(s.c[1], p.c[1])` while for the mixed parameter phase: `s.m = a_T(s.m, p.m)`).
5. `DWPM Hessian.ipynb` [temporary file] contains some of the analytical derivatives which will be written into `nComp.py` functions used to calculate the Hessian soon.
6. `csvDict.py` contains some functions used in the current data handling scheme using the .csv files stored in `\Data`. To be replaced soon.
7. `tgo.py` contains the topographical global optimisation function used to find 
the global minima of many sub-problems in `nComp.py`.

## Config files
The simulation inputs are read from a file called `config.cfg` in the same 
format as `config example.cfg` and needs to placed in the directory before 
running the `main.py` file. The `[dir]` needs to contain a directory to a Data 
folder containing the `.csv` files used in the project. Any `.csv` file with 
the required headings as read in `data_handling.py` can be read.

## Running unittests
The `*_tests.py` files contain unittests. To run unittests covering all 
the code run `test_all.py`.

Specific test suites in each file can be easily run from the command line as 
example:

`$ python2 -m unittest -v tgo_tests.TestTgoFuncs`

Specific function tests inside a unittest can be run for example:

`$ python2 -m unittest -v tgo_tests.TestTgoSubFuncs.test_t1`

NOTE: Only tgo_tests.py are fully working yet, 

