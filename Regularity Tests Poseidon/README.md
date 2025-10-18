# Regularity Tests For Poseidon

### Files:
1) `generate_params_poseidon.sage` is obtained from the [code](https://extgit.isec.tugraz.at/krypto/hadeshash/-/tree/master/code?ref_type=heads) provided by the authors. This script is modified slightly to get desired formatting in the `poseidon_parameters.txt` file
2) `run_regularity_test.py` reads the parameter and calls MAGMA subprocess to compute the Hilbert series, then compares the series to check if they match
3) `regularity.py` is the core computation. It contains the Poseidon round functions. Here we construct the polynomial system, then compute the Hilbert series of the ideal and the expected Hilbert series, both of which are written in files for the python script to read compare.

### Running the experiment

First Run `generate_params_poseidon.sage` with the appropriate command line parameters (These can be found inside the comments in sage file). Four instances are included: 64 bit prime, 254 bit prime with 3 branches, 254 bit prime with 5 branches, 255 bit prime with 3 branches, 255 bit prime with 5 branches 
Then run the python file, specifying the parametrs in the command line seperated by underscores ```--parameters <Prime size in bits>_<RF/2>_<RP>_alpha_t_k``` 

for example
``` 
run_regularity_test.py --parameters 64_2_20_3_24_2
```

Valid argunents for the python script are:
-  Prime size in bits: Can be 64, 254, 255. 
-  $R_F/2$: (Number of Full Rounds/2): Any integer > 1  
- $R_P$: (Number of Partial rounds) Any integer > 1
- alpha : (S-Box exponent) : 3 or 5 
- t : (Number of branches) 3, 5, 24 
- $k$ : (Specifies CICO-$k$ instance) $1 \leq k \leq t$

**Note** Prime size in bits and Number of branches must the same as the parameters generated with the sage script in the first step. Other parameters can virtually take any values but this affects the run time so only some are realistic.