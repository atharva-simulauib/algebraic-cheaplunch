import hashlib
import itertools
from time import time
import sys

p = 2**64 - 2**32 + 1
p = next_prime(p) # for more possible alpha
F = FiniteField(p)

class NoSolutions(Exception):
    pass

def matrix(b):
    matrix_row = [7, 23, 8, 26, 13, 10, 9, 7, 6, 22, 21, 8][:b]
    M = Matrix(F, [[matrix_row[(i-j)%b] for j in range(b)] for i in range(b)])
    return M

def mds_matrix(b):
    x = [F.random_element() for i in range(b)]
    y = [F.random_element() for i in range(b)]
    return Matrix(F, [[1/(xx-yy) for yy in y] for xx in x])

def compute_round_constants(p, m, c, lam, N):
    # Step 1: Format the input string
    input_str = f"RPO({p},{m},{c},{lam})"
    input_bytes = input_str.encode('ascii')
    
    # Step 2: Calculate total number of bytes needed
    num_constants = m * N
    total_bytes = num_constants * 9  # 9 bytes per constant
    
    # Step 3: Expand input string using SHAKE256
    shake = hashlib.shake_256()
    shake.update(input_bytes)
    random_bytes = shake.digest(total_bytes)
    
    # Step 4: Split into 9-byte chunks, interpret as little-endian integers, reduce mod p
    constants = []
    for i in range(num_constants):
        chunk = random_bytes[i*9:(i+1)*9]
        value = int.from_bytes(chunk, byteorder='little')
        constants.append(value % p)
    
    return constants

def degree_from_leading_monomials(leading_mons):
    """
    Compute the degree of a zero-dimensional ideal from the list of leading monomials.

    INPUT:
        leading_mons -- list of monomials (in the polynomial ring) which are the leading monomials of a Groebner basis of an ideal of dimension 1

    OUTPUT:
        The number of monomials not divisible by any of the given leading monomials (degree of the ideal).
    """
    R = leading_mons[0].parent()
    n = R.ngens()
    
    # Convert leading monomials to exponent tuples
    lm_exps = [m.exponents()[0] for m in leading_mons]
    
    # Compute bounding box for exponents based on max exponents in LMs
    max_exps = [max(lm[i] for lm in lm_exps) for i in range(n)]
    
    degree_count = 0

    # Check if a monomial exponent vector is divisible by any leading monomial exponent vector
    def divisible_by_any(exp):
        for lm_exp in lm_exps:
            # lm divides monomial if for all vars lm_exp_i <= exp_i 
            if all(exp[i] >= lm_exp[i] for i in range(n)):
                return True
        return False

    # Recursively enumerate exponents lex order <= max_exps, count if not divisible
    def count_monomials(pos=0, current_exp=[0]*n):
        nonlocal degree_count
        if pos == n:
            if not divisible_by_any(tuple(current_exp)):
                degree_count += 1
            return
        for e in range(max_exps[pos]):
            current_exp[pos] = e
            count_monomials(pos+1, current_exp)
        # Also consider exponent = max_exps[pos]
        current_exp[pos] = max_exps[pos]
        count_monomials(pos+1, current_exp)

    count_monomials()
    return degree_count

def grevlex_gb_and_di(eqs, magma=False):
    """
        Computes a gb of ideal(eqs) using magma and,
        if the ideal has dimension 0,
        return D_I (computed efficiently, don't need lex),
        else return D_I = None
    """
    I = ideal(eqs)
    if magma:
        G = I.groebner_basis(algorithm = "magma")
    else:
        G = I.groebner_basis()
    if ideal(G).dimension() == 0:
        d_i = degree_from_leading_monomials([g.lm() for g in G])
    else:
        d_i = None
    return G, d_i

def magma_variety(eqs):
    # this workaround is because sage won't factor polynomials for big characteristics
    I = ideal(eqs)
    Ilex = I.change_ring(I.ring().change_ring(order='lex'))
    Glex = Ilex.groebner_basis(algorithm="magma")
    return ideal(Glex).variety()



def power7_map_on_fp3(alpha):
    # Define extension field F_{p^3} using modulus X^3 - X - 1
    # (note: the xhash paper doesn't say which modulus they used,
    # i had to reverse engineer it...)
    SR.<s1,s2,s3> = F[]
    FR.<z> = PolynomialRing(SR)
    K.<a> = SR.extension(z^3 - z - 1)
    elem = s1 + s2*a + s3*a^2
    elem_p7 = elem^alpha
    return elem_p7

def round_skip_compute_vector(M, pos, posS, alpha, c):
    # given a matrix M, k (for CICO-k) and a vector of positions (a support for the output),
    # and positions for skipped sboxes on the second step,
    # compute the normalized (1 on first nonzero coef) vector that has this support
    # and verifies the input constraints of the CICO (no round constant additions).
    assert(sum(pos) == c+1 + sum(posS))
    b = M.nrows()
    assert(len(pos) == b)
    assert(len(posS) == b)
    temp_ring = PolynomialRing(M.base_ring(), c + sum(posS), names="a")
    vec = zero_vector(temp_ring, b)
    coefs = [1] + list(temp_ring.gens()) # the 1 means vec is normalized
    j = 0
    for i in range(b):
        if pos[i] == 1:
            vec[i] = coefs[j]
            j += 1
    M_inv_ring = M.inverse().change_ring(temp_ring)
    vec_before_M = M_inv_ring * vec
    eqs = list(vec_before_M[b-c:b])
    M_ring = M.change_ring(temp_ring)
    vec_after_S_and_M = M * vector(temp_ring, [v^alpha for v in vec])
    eqs.extend(list(v for i,v in enumerate(vec_after_S_and_M) if posS[i] == 1)) 
    sols = magma_variety(eqs)
    if len(sols) == 0:
        raise NoSolutions("no solutions, try again!")
    sol = sols[0]
    return vector(M.base_ring(), [v.substitute(sol) for v in vec])

def round_skip_compute_vector_constants(M, pos, c, rc):
    # given a matrix M, k (for CICO-k) and a vector of positions (a support),
    # compute the vector that verifies the input constraints of the CICO 
    # (with round constant additions).
    assert(sum(pos) == c)
    b = M.nrows()
    assert(len(pos) == b)
    temp_ring = PolynomialRing(M.base_ring(), c, names="a")
    vec = zero_vector(temp_ring, b)
    coefs = list(temp_ring.gens()) 
    j = 0
    for i in range(b):
        if pos[i] == 1:
            vec[i] = coefs[j]
            j += 1
    M_inv_ring = M.inverse().change_ring(temp_ring)
    vec_before_M = M_inv_ring * vec - vector(temp_ring, list(rc[:b]))
    eqs = list(vec_before_M[b-c:b])
    sol = ideal(eqs).variety()[0]
    return vector(M.base_ring(), [v.substitute(sol) for v in vec])



class xhash8:
    def __init__(self, b, alpha, num_steps, num_cico_vars, num_inv_alpha=None):
        # the last step is always the MC step
        self.b = b
        if num_inv_alpha == None:
            self.num_inv_alpha = (2*self.b+1)//3 #xhash8
        else:
            self.num_inv_alpha = num_inv_alpha
        self.alpha = alpha
        self.alpha_inv = inverse_mod(alpha, p-1)
        self.N = num_steps
        self.rc = compute_round_constants(p, b, 1, 128, 3*self.N+1)
        self.k = num_cico_vars
        #self.num_partial_vars = self.N * (2*self.b+1)//3
        self.num_partial_vars = self.N * self.num_inv_alpha
        names = [f"y{i+1}" for i in range(self.num_partial_vars)] + [f"x{i+1}" for i in range(self.k)]
        #names = [f"x{i+1}" for i in range(self.k)] + [f"y{i+1}" for i in range(self.num_partial_vars)]

        weights = [1]*(self.k) # x variables
        # y variables weights
        weights.extend([(self.alpha)^(i//(self.num_inv_alpha)) for i in range((self.N) * self.num_inv_alpha)])
        assert len(weights) == len(names)
        weights = [int(w) for w in weights]
        self.can_use_magma = (self.N == 1) # gets overridden if we skip rounds
        if not self.can_use_magma:
            self.R = PolynomialRing(F, self.k + self.num_partial_vars, names=names,
                                       order=TermOrder('wdegrevlex', weights))
        else:
            self.R = PolynomialRing(F, self.k + self.num_partial_vars, names=names, order='degrevlex')

        self.x = self.R.gens()[self.num_partial_vars:]
        self.y = self.R.gens()[:self.num_partial_vars]
        #self.x = self.R.gens()[:self.k]
        #self.y = self.R.gens()[self.k:]
        self.state = vector(self.R, list(self.x) + [0]*(b - self.k))
        self.step_num = 0
        self.eq_list = []
        self.M = mds_matrix(self.b)
        self.p0, self.p1, self.p2 = power7_map_on_fp3(alpha)

    def add_rc(self, nr):
        self.state += vector(self.state.base_ring(), self.rc[nr*self.b : (nr+1)*self.b])

    def inverse_alpha_round(self, step):
        var_num = 0
        for i in [j for j in range(self.b) if j%3 != 1]:
            new_var = self.y[step*(2*self.b+1)//3 + var_num]
            var_num += 1
            self.eq_list.append(new_var^self.alpha - self.state[i])
            self.state[i] = new_var

    def inverse_alpha_round_partial(self, step):
        # not the real inverse round function: for testing purposes
        for i in range(self.num_inv_alpha):
            new_var = self.y[step*self.num_inv_alpha + i]
            self.eq_list.append(new_var^self.alpha - self.state[i])
            self.state[i] = new_var

    def inverse_alpha_round_partial_skip(self, skip_y, num_skipped_sboxes, step):
        # not the real inverse round function: for testing purposes
        actual_skipped_sboxes = min(num_skipped_sboxes, self.num_inv_alpha)
        if step == 0:
            for i in range(actual_skipped_sboxes):
                assert self.state[i].degree() == 0
                self.state[i] = self.state[i]^(self.alpha_inv)
            for i in range(actual_skipped_sboxes, self.num_inv_alpha):
                new_var = skip_y[i - actual_skipped_sboxes]
                self.eq_list.append(new_var^self.alpha - self.state[i])
                self.state[i] = new_var
        else:
            for i in range(self.num_inv_alpha):
                new_var = skip_y[step*self.num_inv_alpha + i - actual_skipped_sboxes]
                self.eq_list.append(new_var^self.alpha - self.state[i])
                self.state[i] = new_var

    def alpha_round(self):
        for i in range(self.b):
            self.state[i] = (self.state[i]^(self.alpha)).reduce(self.eq_list)


    def p3_layer(self):
        for i in range(self.b//3):
            x0,x1,x2 = self.state[3*i:3*i+3]
            y0 = self.p0(x0, x1, x2)
            y1 = self.p1(x0, x1, x2)
            y2 = self.p2(x0, x1, x2)
            self.state[3*i:3*i+3] = (y0, y1, y2)
        for i in range(3*(self.b//3), self.b):
            self.state[i] = self.state[i]^(self.alpha)

    def linear_layer(self):
        self.state = self.M.change_ring(self.state.base_ring()) * self.state

    def compute_eqs(self):
        self.add_rc(0)
        self.linear_layer()
        for step in range(self.N):
            self.alpha_round()
            self.linear_layer()
            self.add_rc(3*step+1)
            self.inverse_alpha_round_partial(step)
            self.add_rc(3*step+2)
            self.p3_layer()
            self.linear_layer()
            self.add_rc(3*step+3)
        self.eq_list.extend(self.state[:self.k])
        return self.eq_list

    def compute_eqs_round_skip(self, verbose=True):
        # From Appendix "Case of Poseidon"
        s = floor((self.b - self.k)/(self.k)) - 1
        num_skipping_variables = min(s, self.k)
        num_skipped_sboxes = max(0, s-self.k)
        actual_skipped_sboxes = min(num_skipped_sboxes, self.num_inv_alpha)
        self.can_use_magma = (self.N == 1) and (num_skipping_variables == 0 or num_skipping_variables == self.k)
        print(f"round skips: {num_skipping_variables} skipping variables, {num_skipped_sboxes} skipped sboxes\n")
        # watch out: there may be more "skipped sboxes" than there are sboxes.
        # we still solve as if there were num_skipped_sboxes to preserve dimension-0 set of solutions.
        skip_partial_vars = (self.N - 1) * self.num_inv_alpha + max(0, self.num_inv_alpha - num_skipped_sboxes)
        names_x = [f"x{i+1}" for i in range(self.k)]
        names_y = [f"y{i+1}" for i in range(skip_partial_vars)]
        assert(len(names_x + names_y) == self.k + skip_partial_vars)

        #give higher weights to variables that skip
        weights_x = [self.alpha]*num_skipping_variables + [1]*(self.k - num_skipping_variables) # x variables
        weights_y = [1]*(self.num_inv_alpha - actual_skipped_sboxes) + [(self.alpha)^(i//(self.num_inv_alpha) + 1) for i in range((self.N - 1) * self.num_inv_alpha)]
        assert len(weights_x) + len(weights_y) == len(names_x) + len(names_y)

        # choose ordering such that y^alpha lead the first equations
        # and x vars leads the last ones
        if not self.can_use_magma:
            skip_ring = PolynomialRing(F, self.k + skip_partial_vars, names=names_x + names_y,
                                       order=TermOrder('wdegrevlex', weights_x + weights_y))
            skip_x = skip_ring.gens()[:self.k]
            skip_y = skip_ring.gens()[self.k:]
            self.weights = weights_x + weights_y
        else:
            if num_skipping_variables == 0:
                skip_ring = PolynomialRing(F, self.k + skip_partial_vars, names=names_y + names_x)
                skip_x = skip_ring.gens()[skip_partial_vars:]
                skip_y = skip_ring.gens()[:skip_partial_vars]
            else:
                skip_ring = PolynomialRing(F, self.k + skip_partial_vars, names=names_x + names_y)
                skip_x = skip_ring.gens()[:self.k]
                skip_y = skip_ring.gens()[self.k:]
        self.x = skip_x
        self.y = skip_y
        self.R = skip_ring

        # now, compute the vectors that represent the state
        # before the first sbox, start with constant part
        skip_state = round_skip_compute_vector_constants(self.M,
                                                         [0]*(self.b - self.k) + [1]*self.k,
                                                         self.k, self.rc)
        
        # now, the variables that skip the first step
        hamming_pos = self.k + 1 + num_skipped_sboxes
        posS = [1]*(num_skipped_sboxes) + [0]*(self.b - num_skipped_sboxes) # assumes simplified P'
        for i in range(num_skipping_variables):
            pos = [0]*(i*hamming_pos) + [1]*(hamming_pos) + [0]*(self.b - (i+1)*hamming_pos)
            vec_xi = round_skip_compute_vector(self.M, pos, posS, self.alpha, self.k)
            skip_state += skip_x[i] * vector(skip_ring, vec_xi)
        # now, the variables that don't skip
        last_pos_skipping_variables = num_skipping_variables * hamming_pos
        posS = [0] * self.b
        for i in range(num_skipping_variables, self.k):
            pos = [0] * (self.b - self.k) + [1] * (self.k)
            pos[last_pos_skipping_variables + i - num_skipping_variables] = 1
            vec_xi = round_skip_compute_vector(self.M, pos, posS, self.alpha, self.k)
            skip_state += skip_x[i] * vector(skip_ring, vec_xi)

        M_inv_ring = self.M.inverse().change_ring(skip_ring)
        state_before_M = M_inv_ring * skip_state - vector(skip_ring, list(self.rc[:self.b]))
        skip_state_change_var = zero_vector(skip_ring, self.b)
        for i in range(last_pos_skipping_variables):
            assert(len(skip_state[i].coefficients()) == 1)
            coef = skip_state[i].coefficients()[0]
            mon = skip_state[i].monomials()[0]
            skip_state_change_var[i] = coef^(self.alpha) * mon
        for i in range(last_pos_skipping_variables, self.b):
            skip_state_change_var[i] = skip_state[i]^alpha
        skip_state_chv_after_M = self.M.change_ring(skip_ring) * skip_state_change_var
       
        if verbose:
            print("\ninput state (before change of variable):")
            print(state_before_M)
            print("\nstate after first affine layer (before change of variable):")
            print(skip_state)
            print("\nstate after first nonlinear layer (after change of variable):")
            print(skip_state_change_var)
            print("\nstate after second linear layer (after change of variable):")
            print(skip_state_chv_after_M)

        self.state = skip_state_chv_after_M

        self.add_rc(1)
        self.inverse_alpha_round_partial_skip(skip_y, num_skipped_sboxes, 0)
        self.add_rc(2)
        self.p3_layer()
        self.linear_layer()
        self.add_rc(3)
        for step in range(1,self.N):
            self.alpha_round()
            self.linear_layer()
            self.add_rc(3*step+1)
            self.inverse_alpha_round_partial_skip(skip_y, num_skipped_sboxes, step)
            self.add_rc(3*step+2)
            self.p3_layer()
            self.linear_layer()
            self.add_rc(3*step+3)
        self.eq_list.extend(self.state[:self.k])
        return self.eq_list

    
    def Q_is_regular(self):
        # call after equations are computed. tests if
        # 1. the GB of Q can be computed from the GB of Qtop, and
        # 2. the GB of the whole system is G(Q) \cup P
        Q = self.eq_list[-self.k:]
        P = self.eq_list[:len(self.eq_list)-self.k]
        subR = PolynomialRing(F, self.k, names="x")
        sub_remove_y = dict()
        for i,x in enumerate(self.x):
            sub_remove_y[x] = subR.gens()[i]
        for y in self.y:
            sub_remove_y[y] = 0
        Qtop = ideal([pol.homogeneous_components()[pol.total_degree()].substitute(sub_remove_y) for pol in Q])
        if self.can_use_magma:
            Gtop = Qtop.groebner_basis(algorithm="magma")
            G = ideal(Q).groebner_basis(algorithm="magma")
            G_whole = ideal(self.eq_list).groebner_basis(algorithm="magma")
        else:
            Gtop = Qtop.groebner_basis()
            G = ideal(Q).groebner_basis()
            G_whole = ideal(self.eq_list).groebner_basis()
        fij = [[f.substitute({subR.gens()[i]: self.x[i] for i in range(self.k)}) for f in g.lift(Qtop)] for g in Gtop]
        G_from_fij = [sum(f[i] * Q[i] for i in range(self.k)) for f in fij]
        can_compute_G_from_Qtop = set([g.lm() for g in G_from_fij]) == set([g.lm() for g in G])
        can_compute_GB_whole = set([g.lm() for g in G + P]) == set([g.lm() for g in G_whole])

        return can_compute_G_from_Qtop, can_compute_GB_whole

# to test round skips:
# b = 8, c = 2 for s = c
# b = 6, c = 2 for s < c
# b = 10, c = 2 for s > c
num_steps = 1

def compute_data(b, alpha, num_cico_vars, num_inv_alpha, skip_rounds=False):
    start = time()
    while True:
        try:
            model = xhash8(b, alpha, num_steps, num_cico_vars, num_inv_alpha)
            if skip_rounds:
                eqs = model.compute_eqs_round_skip(verbose=False)
            else:
                eqs = model.compute_eqs()
            break
        except NoSolutions:
            print("no solutions, trying another...")
    I = ideal(eqs)
    G, D_I = grevlex_gb_and_di(eqs, magma=model.can_use_magma)
    deg_reg = max(g.degree() for g in G)
    deg_reg_bound = sum(eq.degree() - 1 for eq in model.eq_list[-model.k:]) + 1
    dibound = prod(eq.total_degree() for eq in model.eq_list)
    can_compute_G_from_Qtop, can_compute_GB_whole = model.Q_is_regular()
    end = time()
    # print(f"took {floor(end - start)}s")
    return D_I, dibound, deg_reg, deg_reg_bound, can_compute_G_from_Qtop, can_compute_GB_whole

def print_line(b, alpha, num_cico_vars, num_inv_alpha):
    di, dib, dr, drb, okbool, _  = compute_data(b, alpha, num_cico_vars, num_inv_alpha)
    if okbool:
        mark = "\\cmark"
    else:
        mark = "\\xmark"
    print(f"({b},{num_cico_vars},{num_inv_alpha}) & {dr} & ({drb}) & {di} & ({dib}) & {mark} \\\\") 

experiment_bound = 2^10
for alpha in [3,5]:
    for b in range(2,12):
        for num_cico_vars in range(2,b+1):
            for num_inv_alpha in range(1,b):
                dib = alpha^(2*num_cico_vars + num_inv_alpha)
                if dib < experiment_bound:
                    # print(f"b = {b}, num_inv_alpha = {num_inv_alpha}")
                    print_line(b, alpha, num_cico_vars, num_inv_alpha)



#b = 8
#alpha = 3
#num_steps = 1
#num_cico_vars = 2
#num_inv_alpha = 1
#print("computing equations...")
#while True:
#    try:
#        model = xhash8(b, alpha, num_steps, num_cico_vars, num_inv_alpha)
#        #eqs = model.compute_eqs_round_skip(verbose=True)
#        eqs = model.compute_eqs()
#        break
#    except NoSolutions:
#        print("no solutions, trying another...")
#print(f"USING MAGMA: {model.can_use_magma}")
#print("computing groebner basis...")
#I = ideal(eqs)
#G, D_I = grevlex_gb_and_di(eqs, magma=model.can_use_magma)
#print("degree of regularity:")
#print(max(g.degree() for g in G))
#print("wdreg bound on Qtop:")
#print(sum(eq.degree() - 1 for eq in model.eq_list[-model.k:]) + 1)
#dibound = prod(eq.total_degree() for eq in model.eq_list)
#print(dibound)
#print("degree of ideal:\n", D_I, f"(upper bound: {dibound})")
##Ilex = I.change_ring(I.ring().change_ring(order='lex'))
##Glex = Ilex.groebner_basis(algorithm="magma")
##print("ideal degree:", Glex[-1].degree(), " = ", factor(Glex[-1].degree()))
#can_compute_G_from_Qtop, can_compute_GB_whole = model.Q_is_regular()
#print("can compute G(Q) from Qtop:", can_compute_G_from_Qtop)
#print("can compute G(P,Q) as P cup G(Q):", can_compute_G_from_Qtop)


#R = model.R
#x = model.x
#y = model.y
#subring = PolynomialRing(F, num_cico_vars, names=[f"x{i+1}" for i in range(num_cico_vars)])
#sr = subring.gens()
#sub_dic = dict()
#for xi, sri in zip(x, sr):
#    sub_dic[xi] = sri
#for yi in y:
#    sub_dic[yi] = 0

#top_homo_eqs = [pol.homogeneous_components()[pol.total_degree()] for pol in eqs[-num_cico_vars:]]
#top_homo_ideal = ideal(top_homo_eqs)
#GH = top_homo_ideal.groebner_basis()
#print("The leading monomials of the GB (top homogeneous parts of k equations) are:")
#print([g.lm() for g in GH])
#
#y = model.y
#eqs_only_x = [sum(pol.monomial_coefficient(m)*m for m in pol.monomials() if all(m.gcd(yi)==1 for yi in y)) for pol in eqs[-num_cico_vars:]]
#ideal_only_x = ideal(eqs_only_x)
#GoX = ideal_only_x.groebner_basis()
#print("The leading monomials of the GB (only x monomials) are:")
#print([g.lm() for g in GoX])
