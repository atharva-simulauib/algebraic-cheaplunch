# Regularity tests for Neptune

### Files:

1) `generate_mds.sage` generates random MDS matrices for the external rounds following the specific design as in https://doi.org/10.46586/tosc.v2022.i3.20-72
2) `run_regularity_test.py` accepts one argument --parametrs in which is a string in the format  "<Prime size in bits> <RF1> <RP> <RF2> t k". It calls sage subprocess to generate random MDS matrix for this particular instancee,then  calls MAGMA subprocess to generate the polynomial system for CICO-$k$ modelling and compute the Hilbert series, then compares the series to check if they match
3) `regularity.py` is the core computation. It contains the Neptune round functions. Here we construct the polynomial system, then compute the Hilbert series of the ideal and the expected Hilbert series, both of which are written in files for the python script to read compare.

### Running the experiment 

Run the python file `run_regularity_test.py` with argument --parametrs in which is a string in the format  "<Prime size in bits> <RF1> <RP> <RF2> t k". 

For example 
```
python3 run_regularity_tests.py --parameters "0xffffffffffffffc5 2 6 1 3 2"
```
This will run the experiment for CICO-2 problem for Neptune in 64.bit prime field  with 6 branches and $RF1=2. RP=6, RF2=1$

Other primes tested:
BLS12 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

BN254 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
