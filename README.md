[![Code Climate](https://codeclimate.com/github/Stefan-Endres/DWPM-Mixture-Model/badges/gpa.svg)](https://codeclimate.com/github/Stefan-Endres/DWPM-Mixture-Model)
[![Test Coverage](https://codeclimate.com/repos/56b4f4a0983f2b007d00135e/badges/
71a1b441b1ae6f925618/coverage.svg)](https://codeclimate.com/repos/
56b4f4a0983f2b007d00135e/coverage)
[![Issue 
Count](https://codeclimate.com/repos/56b4f4a0983f2b007d00135e/badges/
71a1b441b1ae6f925618/issue_count.svg)](https://codeclimate.com/repos/
56b4f4a0983f2b007d00135e/feed)


# DWPM-Mixture-Model
Phase seperation calculation using the DWPM mixture rule.

## Main function call examples.
The functions can be used by calling the `main.py` script with optionals 
passed from the command line, for example initializing the state and parameter 
classes `s` `p` for acetone water phenol with two valid root phases and 
then enter the Python shell:

`$ python -i main.py -c acetone benzene water -p x y`
```python
>>> s
<nComp.state instance at 0x7f28190cf368>
>>> p
<data_handling.MixParameters instance at 0x7f28190d01b8>

```

[To be implemented soon:]
To optmise new paramters and save the results add the `-optim` `-save` options

`$ python main.py -c acetone benzene water -p x y -optim -save` 

To calculate all phase seperations and equilibrium points specifify the 
pressure and temperature only

`$ python main.py -c acetone water -p x y -P 101e3 -T 298.15`

To calculate equilibrium at a specific feed point specifify the pressure, 
temperature and a feed composition

`$ python main.py -c acetone water -p x y -P 101e3 -T 298.15 -z 0.5 0.4`

Several plots for binary and ternary systems are available.
The plot optional `-pltg` can be added to plot the Gibbs surface and any 
equilibrium tangent planes at a P, T point. 

`$ python main.py -c acetone water -p x y -P 101e3 -T 298.15 -pltg `

The `-pltiso` will plot either isotherms or isobars with the current optimised 
parameters at the data $P$, $T$, $\bar{x}$ data points to compare the fit

`$ python main.py -c acetone water -p x y -pltiso`

Specifying `-pltiso` with a temperature or pressure point will plot an isotherm 
or isotherm at that point

`$ python main.py -c acetone water -p x y -T 298.15 -pltiso`

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

