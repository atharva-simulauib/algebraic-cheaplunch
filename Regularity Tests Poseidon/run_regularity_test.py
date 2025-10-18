import subprocess
import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument('--parameters', type=str, help="Poseidon parametrs fot the regularity test in the form of <Prime size in bits>_<RF/2>_<RP>_alpha_t_k")
args = parser.parse_args()

input_parameters = [int(j) for j in args.parameters.split('_')]
output_file= open( f"{args.parameters}_output.txt","w")
sys.stdout = output_file


print("\n--> Reading parametrs from poseidon_parametrs.txt\n", flush=True)
with open("poseidon_parameters.txt","r") as params_file :
    params_dict = {key.strip(): value.strip() for line in params_file.readlines() for key, value in [line.split(":")]}

print("Suggested parameters: ",params_dict['Params'], flush=True)


magma_input = """
PRIME := {};
gf := FiniteField(PRIME);
RF := {}; // Number of full rounds / 2
RP := {}; // Number of middle partial rounds
alpha := {};        // S-Box exponent
t := {};            // Number of branches
k := {};  // CICO instance (k < t)
MDS := Matrix(gf, t, t, [ StringToInteger(entry[3..#entry],16): entry in &cat {}]);
arg_string := "{}";
""".format(params_dict['Prime number'], 
 input_parameters[1], input_parameters[2], input_parameters[3], input_parameters[4], input_parameters[5],
 params_dict['MDS matrix'].replace("'",'"'), args.parameters )
with open("input.magma", "w") as f:
    f.write(magma_input)

result = subprocess.call(["magma","regularity.magma"], stdout=output_file)
print("\n--> Finished magma computation\n", flush=True)



print("\n--> Comparing Hilbert series\n", flush=True)
expectedHS =open(f"{args.parameters}_HS_expected.txt", "r").read().replace(' ', '').replace('\n','').split('+')[::-1]
systemHS = open(f"{args.parameters}_HS_system.txt", "r").read().replace(' ', '').replace('\n','').split("+")

i = 0
SeriesMatch = True
while i<min(len(expectedHS), len(systemHS)):
    if expectedHS[i+1] != systemHS[i]:
        print(expectedHS[i], systemHS[i])
        SeriesMatch = False
    i += 1
if SeriesMatch: print("\n\n--> Series Match!\n", flush=True)
    