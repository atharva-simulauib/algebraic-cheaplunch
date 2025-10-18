"""
    Reference Implementatopm of ArionHash.
"""

class ArionHash:
    
    def __init__(self,
                 field=GF(0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001), # BLS12
                 branches=3,
                 rounds=6,
                 capacity=1,
                 d_1=None,
                 d_2=None,
                 constants_g=None,
                 constants_h=None,
                 constants_aff=None,
                 initial_value=None):
        self.field = field
        if not self.field.is_prime_field():
            raise Exception("Only prime fields are implemented. " + str(self.field) + " is not a prime field.")
        self.branches = branches
        self.rounds = rounds
        self.capacity = capacity
        if self.capacity < 1 or self.capacity > self.branches:
            raise Exception("Capacity must be greater than zero and smaller than the branch number.")
        self.rate = self.branches - self.capacity
        
        self.d_1 = d_1
        if self.d_1 is None:
            self.d_1 = 2
            while gcd(self.d_1, self.field.characteristic() - 1) != 1:
                self.d_1 += 1
        else:
            if gcd(self.d_1, self.field.characteristic() - 1) != 1:
                raise Exception("d_1 = " + str(self.d_1) + " does not induce a power permutation over " + str(self.field) + ".")
        self.e_1 = xgcd(self.d_1, self.field.characteristic() - 1)[1] % (self.field.characteristic() - 1)
        
        self.d_2 = d_2
        if self.d_2 is None:
            found = False
            for d in [257, 161, 129, 125, 123, 121]:
                if gcd(d, self.field.characteristic() - 1) == 1:
                    self.d_2 = d
                    found = True
                    break
            if not found:
                raise Exception("Could not find power permutation d_2 with exponent in " + str([257, 161, 129, 125, 123, 121])  + " over " + str(self.field) + ".")
        else:
            if gcd(self.d_2, self.field.characteristic() - 1) != 1:
                raise Exception("d_2 = " + str(self.d_2) + " does not induce a power permutation over " + str(self.field) + ".")
        self.e_2 = xgcd(self.d_2, self.field.characteristic() - 1)[1] % (self.field.characteristic() - 1)
        
        self.constants_g = constants_g
        if self.constants_g is None:
            self.constants_g = []
            for r in range(0, self.rounds):
                for k in range(0, self.branches - 1):
                    c_1 = self.field.random_element()
                    c_2 = self.field.random_element()
                    while legendre_symbol(c_1**2 - 4 * c_2, self.field.characteristic()) != -1:
                        c_1 = self.field.random_element()
                        c_2 = self.field.random_element()
                    self.constants_g.append([c_1, c_2])
        else:
            if len(flatten(self.constants_g)) != self.rounds * (self.branches - 1) * 2:
                raise Exception("Number of round constants for the polynomials g_i does not match branches minus one times the number of rounds times two.")
            for pair in self.constants_g:
                if legendre_symbol(pair[0]**2 - 4 * pair[1], self.field.characteristic()) != -1:
                    raise Exception("The constant pair (" + str(pair[0]) + ", " + str(pair[1]) + ") for g_i has zeros over " + str(self.field) + ".")
        
        self.constants_h = constants_h
        if self.constants_h is None:
            self.constants_h = []
            for i in range(0, self.rounds * (self.branches - 1)):
                self.constants_h.append(self.field.random_element())
        else:
            if len(self.constants_h) != self.rounds * (self.branches - 1):
                raise Exception("Number of round constants for the polynomials h_i does not match branches minus one times the number of rounds.")
        
        self.constants_aff = constants_aff
        if self.constants_aff is None:
            self.constants_aff = []
            for r in range(0, self.rounds):
                self.constants_aff.append([self.field.random_element() for i in range(0, self.branches)])
        else:
            if len(flatten(self.constants_aff))!= self.branches * self.rounds:
                raise Exception("Number of affine round constants does not match branches times the number of rounds.")
        
        if matrix_type == 1:
            self.matrix = matrix.circulant(vector(self.field, range(1, self.branches + 1)))
        elif matrix_type == 2:
            self.matrix = matrix(self.field, self.branches * [list(range(1, self.branches + 1))])
            for i in range(1, self.branches):
                self.matrix[i, i - 1] += 1
                self.matrix[i, i] -= 1
        else:
            raise Exception("Matrix type not implemented.")
        
        self.initial_value = initial_value
        if self.initial_value is None:
            self.initial_value = self.capacity * [0]
        else:
            if len(initial_value) != self.capacity:
                raise Exception("Length of initial value does not match sponge capacity.")
        
        print("ArionHash parameters")
        print("Prime field: " + str(self.field.characteristic()))
        print("Branches: " + str(self.branches))
        print("Rounds: " + str(self.rounds))
        print("Sponge capacity: " + str(self.capacity))
        print("Exponent d_1: " + str(self.d_1))
        print("Exponent d_2: " + str(self.d_2))
        print("Matrix:\n" + str(self.matrix))
        print("Constants for the g_i's: " + str(self.constants_g))
        print("Constants for the h_i's: " + str(self.constants_h))
        print("Affine constants: " + str(self.constants_aff))
        print("Initial value: " + str(self.initial_value))
    
    def evaluate_g_and_h(self, x_in, constants_g, constant_h):
        # Constant term
        out_g = constants_g[1]
        out_h = self.field(0)
        # Linear term
        out_g += x_in * constants_g[0]
        out_h += x_in * constant_h
        # Quadratic term
        x_temp = x_in**2
        out_g += x_temp
        out_h += x_temp
        return out_g, out_h
    
    def GTDS(self, v_in, constants_g, constants_h):
        v_out = copy(v_in)
        v_out[self.branches - 1] = v_out[self.branches - 1]**self.e_2
        sigma = v_in[self.branches - 1] + v_out[self.branches - 1]
        for i in range(self.branches - 2, -1, -1):
            v_out[i] = v_out[i]**self.d_1
            g_i, h_i = self.evaluate_g_and_h(sigma, constants_g[i], constants_h[i])
            v_out[i] *= g_i
            v_out[i] += h_i
            sigma += v_in[i] + v_out[i]
        return v_out
    
    def affine_layer(self, v_in, round_constants):
        return self.matrix * v_in + round_constants
    
    def round_function(self, v_in, constants_g, constants_h, constants_aff):
        v_out = self.GTDS(v_in, constants_g, constants_h)
        v_out = self.affine_layer(v_out, constants_aff)
        return v_out
    
    def hash(self, plain):
        if not type(plain) is list:
            plain = [plain]
        # Padding of plain
        m = len(plain)
        initial_value = copy(self.initial_value)
        if len(plain) % self.rate != 0:
            while len(plain) % self.rate != 0:
                plain.insert(0, 0)
            initial_value[0] = m
        initial_value = [self.field(el) for el in initial_value]
        state = vector(self.field, self.rate * [0] + initial_value)
        # Hashing of plain
        for i in range(0, len(plain), self.rate):
            for j in range(0, self.rate):
                state[j] += plain[i:i + self.rate][j]
            state = self.matrix * state
            for r in range(0, self.rounds):
                state = self.round_function(state,
                                            self.constants_g[r * (self.branches - 1):(r + 1) * (self.branches - 1)],
                                            self.constants_h[r * (self.branches - 1):(r + 1) * (self.branches - 1)],
                                            vector(self.field, self.constants_aff[r]))
        return state[0]