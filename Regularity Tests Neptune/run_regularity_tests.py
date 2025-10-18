import subprocess
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--parameters', type=str, help="Neptune parametrs fot the regularity test in the form of <Prime size in bits> <RF1> <RP> <RF2> t k")
args = parser.parse_args()
arg_string = (args.parameters).replace(" ", "_")
print(arg_string)
output_file= open( f"{arg_string}_output.txt","w")
sys.stdout = output_file

subprocess.call(["sage", "generate_mds.sage", args.parameters])

result = subprocess.call(["magma","regularity.magma"], stdout=output_file)
print("\n--> Finished magma computation\n", flush=True)

print("\n--> Comparing Hilbert series\n", flush=True)
expectedHS =open(f"{arg_string}_HS_expected.txt", "r").read().replace(' ', '').replace('\n','').split('+')[::-1]
systemHS = open(f"{arg_string}_HS_system.txt", "r").read().replace(' ', '').replace('\n','').split("+")

i = 0
SeriesMatch = True
while i<min(len(expectedHS), len(systemHS)):
    if expectedHS[i+1] != systemHS[i]:
        print(expectedHS[i], systemHS[i])
        SeriesMatch = False
    i += 1
if SeriesMatch: print("\n\n--> Series Match!\n", flush=True)
    
