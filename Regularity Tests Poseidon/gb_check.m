// Poseidon regularity test in Magma
SetMemoryLimit(40000000000);
// load the generated parameters
load "input.magma";

// redefine small RF and RP 
RF := 2;
RP := 2;
printf "Experiment parameters: ";
printf "Field with modulus %o\n", PRIME;
printf "RF=%o, RP=%o, alpha=%o, t=%o, CICO-%o instance\n", RF, RP, alpha, t, k;
printf "MDS dimension: %o\n", NumberOfRows(MDS);
SET_CONSTANT_ZERO := false;                     // skip constant addition, for getting only highest degree part

// ------------------------------------  Poseidon Round Function ------------------------------------
// Choose distinct x_i and y_j
coeff_x := [Random(gf) : i in [1..t]];
coeff_y := [Random(gf) : i in [1..t]];
// Define the MDS matrix with entries Mi,j = 1 / (x_i + y_j)
// Uncomment the following line to redefine  MDS as a random matrix (instead of using from Poseidon parameters)
//MDS := Matrix(gf, t, t, [1 / (coeff_x[i] + coeff_y[j]) : i in [1..t], j in [1..t]]); 

// Poseidon round function components
function AddC(in_state)
    if SET_CONSTANT_ZERO eq true then
        return in_state;
    end if;
    return [in_state[i] + Random(gf) : i in [1..#in_state]];
end function;

function MatMul(in_state)
    res := [ &+[ in_state[i] * MDS[i,j] : i in [1..t] ] : j in [1..t] ];
    return res;
end function;

function Sbox(in_state)
    return [in_state[i]^alpha : i in [1..#in_state]];
end function;

function Partial_Sbox(in_state)
    res := in_state;
    res[1] := in_state[1]^alpha;
    return res;
end function;

// Custom functions for debugging and printing
procedure PrintPolys(polyList)
    for i in [1..#polyList] do
        printf "Polynomial %o: with degrees %o\n", i, [Degree(polyList[i], j) : j in [1..k+RP]];
    end for;
    printf "\n";
end procedure; 

// ------------------------------------  Polynomial ring and system setup ------------------------------------
// Ring = Fp[y_1,y_2,...,y_RP, x_1,..., x_k]
weight_vector := [alpha^(RF-1) : i in [1..RP]] cat [1 : i in [1..k]]; 
R := PolynomialRing(gf, weight_vector);
printf "\n";
// Define variables
var_names := ["y" cat IntegerToString(j) : j in [1..RP]] cat
             ["x" cat IntegerToString(i) : i in [1..k]];
AssignNames(~R, var_names);

vars := [R.i : i in [1..k+RP]];     // Extract variables from the ring
R; 

print "Monomial order check: y1 > y2 > ... y_RP > x_1^(alpha^(RF-1)) > x_k^(alpha^(RF-1)) " ;
printf "%o\n",[R.(j) gt R.(j+1) : j in [1..RP]] cat [ R.(RP) gt R.(RP+1)^((alpha^(RF-1))) ] cat [R.(j)^((alpha^(RF-1))) gt R.(j+1)^((alpha^(RF-1))) : j in [RP+1..RP+k-1]];

// ------------------------------------ Modelling ------------------------------------
state := vars[RP+1..RP+k] cat [R!0 : i in [1..t-k]];                        // Initial state containing [x1,..,xk,0,...0]   
printf "Initial state: %o\n", state;
polySys := [];                                                              // Polynomial system where we keep track of new relations

// First full rounds
for j in [1..RF] do
    printf "First full round %o\n", j;                                                         
    state := Sbox(MatMul(AddC(state)));
end for;
state := MatMul(state);
printf "After first full rounds state is:\n";
PrintPolys(state);

// Partial rounds with reduction
for j in [1..RP] do
    printf "Partial round %o\n", j;                                                         
    Append(~polySys, vars[j] - state[1]);                                   // Define yi - state[1] = 0 for each partial round, add the relation to polySys
    state := [vars[j]] cat state[2..t];                                     // Update the first branch with the new variable
    state := MatMul(Partial_Sbox(AddC(state)));

    // reducing the state
    S := quo<R | Ideal(polySys)>;                                           // R / <P1,..,Pi-1>
    //print "%o\n", S;
    for i in [1..#state] do
        state[i] := ChangeRing(state[i],S);
    end for;

end for;

// Final full rounds with reduction
for j in [1..RF] do
    printf "Last Full round %o\n", j;                                                         
    state := Sbox(MatMul(AddC(state)));
    for i in [1..#state] do                               
        state[i] := ChangeRing(state[i],S);
    end for;
end for;

printf "Final state is:\n";
PrintPolys(state);

// Add equations Q1, Q2, ..., Qk = 0 to polySys  (final CICO-k constraints)
for i in [1..k] do
    Append(~polySys, state[i]);
end for;

// Final system contains P1,Q1,Q2,...,QK (in reduced form)
printf "System contains:\n";
for i in [1..#polySys] do
    printf "Polynomial %o: with degrees %o, LT: %o\n", i, [Degree(polySys[i], j) : j in [1..k+RP]], LeadingTerm(polySys[i]);
end for;


freeSys := [polySys[i] : i in [2..RP]];
cheapSys := [polySys[1]] cat [polySys[i] : i in [RP+1..#polySys] ];      
homSys :=  [ HomogeneousComponent(cheapSys[i], WeightedDegree(cheapSys[i]))  : i in [1..#cheapSys] ];
J := IdealWithFixedBasis(homSys);
printf "Computing Coordinate MAtrix of Qtop\n";
MatJ:= CoordinateMatrix(J);

printf "GB {Qtop} -> GB {P,Q}\n";
gb := freeSys;
for i in [1..NumberOfRows(MatJ)] do
    gb_element := &+[ MatJ[i,j] * cheapSys[j] : j in [1..#cheapSys] ] ;
    print LeadingMonomial(gb_element);
    Append(~gb, gb_element);
end for;
SetVerbose("Groebner",2);
printf "<LM(I)>\n:%o", LeadingMonomialIdeal(Ideal(polySys));
printf "Is a Groebner Basis %o\n", LeadingMonomialIdeal(Ideal(polySys)) eq LeadingMonomialIdeal(Ideal(gb));