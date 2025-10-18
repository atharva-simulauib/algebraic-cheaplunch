import sys

# Get the input argument
arg = sys.argv[1]
input_parameters = [int(arg.split(' ')[0], 16)] + [int(j) for j in arg.split(' ')[1:]]
arg_string = arg.replace(" ", "_")

PRIME, RF1, RP, RF2, t_dash, k = input_parameters
t = 2*t_dash
WRITE_MAGMA_INPUT = True

M_ = matrix.circulant( random_vector(GF(PRIME), t_dash) )
M__ = matrix.circulant( random_vector(GF(PRIME), t_dash) )

M_ext = matrix(GF(PRIME), t)

for i_ in range(t_dash):
    for j_ in range(t_dash):
        M_ext[2*i_, 2*j_] = M_[i_, j_]
for i__ in range(t_dash):
    for j__ in range(t_dash):
        M_ext[2*i__ + 1, 2*j__ +1] = M__[i__, j__]

print(M_)
print(M__)
print(M_ext)
if WRITE_MAGMA_INPUT :
    matrix_string = "["
    for i in range(0, t):
        matrix_string += str([ entry for entry in M_ext[i] ])
        
        if i < (t-1):
            matrix_string += ","
    matrix_string += "]"
    print(matrix_string)
    
    magma_code = f"""PRIME:={PRIME};
    RF1 := {RF1};
    RP := {RP};
    RF2 := {RF2};
    t := {t};
    k := {k};
    t_dash := {t_dash};
    MDS_ext := Matrix(GF({PRIME}), t, t, [ entry: entry in &cat {matrix_string} ]);
    arg_string := "{arg_string}";
    """
    print(magma_code)
    
    with open("input.magma",'w') as FILE:
        FILE.write(magma_code)
    