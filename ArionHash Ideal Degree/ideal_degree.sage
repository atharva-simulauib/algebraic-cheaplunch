load("ArionHash.sage")

def generate_variables(branches, rate, rounds):
    variables = []
    for i in range(0, rate):
        variables.append("x_in_" + str(i + 1))
    variables.append("z_" + str(1))
    for i in range(0, rounds - 1):
        for j in range(0, branches):
            variables.append("x_" + str(j + 1) + "__" + str(i + 1))
        variables.append("z_" + str(i + 2))
    for i in range(1, branches):
        variables.append("x_out_" + str(i + 1))
    return variables

def generate_weight_vector(branches, cico_constraint, rounds, d_1, d_2):
    W = [1]
    for i in range(1,rounds):
        degree_i = ( 2**(branches-1) * (d_1+1) - d_1 )**i 
        W.append( degree_i/ d_2 + 1 )
    for i in range(cico_constraint):
        W.append(1)
    return W

# ArionHash evaluation round by round
def evaluate_g_and_h(x_in, constants_g, constant_h):
    # Constant term
    out_g = constants_g[1]
    out_h = 0
    # Linear term
    out_g += x_in * constants_g[0]
    out_h += x_in * constant_h
    # Quadratic term
    x_temp = x_in**2
    out_g += x_temp
    out_h += x_temp
    return out_g, out_h

def GTDS(v_in, v_n_std, d_1, constants_g, constants_h):
    branches = len(v_in)
    v_out = copy(v_in)
    sigma = v_n_std + v_out[branches - 1]
    for i in range(branches - 2, -1, -1):
        v_out[i] = v_out[i]**d_1
        g_i, h_i = evaluate_g_and_h(sigma, constants_g[i], constants_h[i])
        v_out[i] *= g_i
        v_out[i] += h_i
        sigma += v_in[i] + v_out[i]
    return v_out

def affine_layer(v_in, matrix, round_constants):
    return matrix * v_in + round_constants


def sysGen(branches, cico_constraint, rounds=1, d_1=3, d_2= 7):
    arion_hash = ArionHash(field=field,
                                branches=branches,
                                rounds=rounds,
                                capacity=cico_constraint,
                                d_1=d_1,
                                d_2=d_2,
                                constants_g=None,
                                constants_h=None,
                                constants_aff=None)

    variables =  ["y" + str(i) for i in range(rounds)] + ["x" + str(i) for i in range(cico_constraint)]
    print("Variables:",variables)
    W = generate_weight_vector(branches, cico_constraint, rounds, d_1, d_2)
    denoms = [_.denominator() for _ in W]
    W = [_*LCM(denoms) for _ in W] 
    print("weight: " +str(W))

    P = PolynomialRing(field, variables, order=TermOrder('wdegrevlex',W))
    print(P)
    print("Monomial order: ", P.term_order() )
    variables = [P(var) for var in variables]

    current_state = variables[-cico_constraint:] + [P(0)]*(arion_hash.branches-cico_constraint)
    print("Initial state: ", current_state)
    counter = 0
    polySys = []

    for r in range(rounds):
        current_state = affine_layer( vector(current_state), arion_hash.matrix, 0)
        tmp = current_state[branches - 1]
        current_state = GTDS(current_state, tmp, arion_hash.d_1, arion_hash.constants_g, arion_hash.constants_h)
        polySys.append(variables[counter]**d_2 - tmp)

        current_state[branches - 1] = variables[counter]

        counter += 1    
        print("Round ", r)
        for branch in current_state:
            print(branch.lm())

        if r == (rounds -1) :
            current_state = affine_layer( vector(current_state), arion_hash.matrix, 0)

    polySys += current_state[:cico_constraint]
    print("System contains: ")
    for poly in polySys:
        print("Polynomial with LM:", poly.lm())

    I = ideal(polySys)
    D_I = None
    if compute_ideal_degree:
        D_I = len(I.normal_basis())
        print( "Ideal degree:", D_I )
        if D_I == 1:
            print(I.normal_basis())
        weights_product = int(1)
        for w in  W:
            weights_product *= w 
        delta = lambda t,j : int( 2**(t-j)*(d_1+1) - d_1 )
        print( "Weighted Bezout bound (naive):", ( ( delta(branches, 1)**cico_constraint ) * d_2**rounds ) / weights_product )
        print("Claimed bound:", ( ( delta(branches, 1)* delta(branches, cico_constraint+1)**(cico_constraint-1)  )* d_2**rounds ) / weights_product )
        print("\n\n")
    #return D_I, polySys, P
    return D_I


PRIME = 15013         # PRIME = 10007, d1 = 3 or PRIME = 15013, d1 = 5
field = GF(PRIME)
matrix_type = 1

compute_ideal_degree = True
write_magma_input = False

# set of parameters to test D_I for (in the format (branches, cico_constraint))
parameters_list = [ (4,2), (5,2), (5,3), (6,2)]
D_I_dict = { (arg1, arg2): sysGen(arg1, arg2, d_1=5, d_2=7) for (arg1, arg2) in parameters_list }
print(D_I_dict)

# To test upper bound holds for the tested parameters: 
def arion_hash_ideal_degree_bound(t, k, r=1, e, alpha):
    degree_j = lambda j: 2**(t - j) * (e + 1) - e
    degree_q = 1
    for j in range(1,k+1):
        degree_q *= degree_j(j)
    bezout_bound = alpha * (degree_q)
    return bezout_bound

params = [#[t, k, e, alpha, D_I]
        # p = 10007
        [4,2,3,7, 623],
        [5,2,3,7, 1855],
        [5,3,3,7, 3619],
        [6,2,3,7, 4991],
        [4,2,3,257, 22873],
        # p = 15013
        [4,2,5,7, 1351],
        [5,2,5,7, 4207],
        [5,3,5,7, 11557],
        [6,2,5,7, 11599],
        [4,2,5,257, 49601]
           ]


print("t", "\t",
      "k", "\t",
      "e", "\t",
      "alpha", "\t",
      "D_I", "\t",
      "Bound", "\t"
      )

for param in params:
    t = param[0]
    k = param[1]
    e = param[2]
    alpha = param[3]
    D_I = param[4]
    bound = arion_hash_ideal_degree_bound(t,k,e,alpha)
    
    print(t, "\t",
          k, "\t",
          e, "\t",
          alpha, "\t",
          D_I, "\t",
          bound
          #"True" if (D_I<bound) else "False" 
         )



#D_I, polySys, P = sysGen(branches=4,cico_constraint=2,rounds=2)



'''
if write_magma_input:
    weights = P.term_order().weights()
    var_names = ",".join(P.variable_names())
    weight_str = "[" + ",".join(str(w) for w in weights) + "]"
    magma_input = f""" P<{var_names}> := PolynomialRing(FiniteField({str(PRIME)}), {weight_str});
    I := Ideal({polySys});
    """ 
    with open("input.magma", "w") as f:
        f.write(magma_input)


'''
