# DWPM-Mixture-Model
Phase seperation calculation using the DWPM mixture rule.

### Basic description of files.
1. ```nComp.py``` contains the core functions used to simulate multicomponent equilibrium problems and optimise parameter models. 
2. ```binary.py``` is superseded by ```nComp.py```. The solution methods for the equilibrium calculation are based on equation based approaches (and might not always be at the true global minima). For now it is only used to validate that the functions in ```nComp.py``` are working by comparing some systems with known parameters or giving an expected output (such as matching isotherms etc.). Once all the functions in ```nComp.py``` working and validated this file will be deleted.
3. ```pure.py``` is used to fit single component Van der Waals EOS temperature dependent activity coefficient parameters to the Adachi-Lu and Soave models. These parameters are used in the multicomponent functions. The single component data is contained in ```\Data\pure``` and every component simulated in multicomponent system requires vapour pressure data or parameters stores in the .csv files in ```\Data\pure```. Unlike  ```binary.py``` these functions will always be used.
4. ```Van_der_Waals.py``` contains the functions the volume root, saturation pressure and state dependent functions which are called in both ```pure.py``` and ```nComp.py``` (note that for example, the function ```a_T(s, p)``` is used in calculating the temperature dependent activity coefficient parameters for both the pure and the mixed phases states (for example to calculated for component 1: ```s.c[1] = a_T(s.c[1], p.c[1])``` while for the mixed parameter phase: ```s.m = a_T(s.m, p.m)```).
5. ```DWPM Hessian.ipynb``` [temporary file] contains some of the analytical derivatives which will be written into ```nComp.py``` functions used to calculate the Hessian soon.
6. ```csvDict.py``` contains some functions used in the current data handling scheme using the .csv files stored in ```\Data```. To be replaced soon.
7. ```tgo.py``` contains the topographical global optimisation function used to find the global minima of many sub-problems in ```nComp.py```.
