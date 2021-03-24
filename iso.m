load "vOW4SIKE/magma/arithmetic.m";
load "vOW4SIKE/magma/step.m";
load "edwards.m";

costS := 7.2; //7.2;
costI := 10.8; //10.8;


// ----------------------------- //
// ---- Isogeny Computation ---- //
// ----------------------------- //


// Optain an optimal strategy for traversing the isogeny tree based on relative costs
// p for l-scalar mult and q for l-isogeny computation.
// [Taken from vOW4SIKE]
function OptStrat(n,p,q)
    S := [[]];
    C := [RealField()!0];
    for i in [2..(n)] do
        newCpqs := [(C[i-b] + C[b] + b*p + (i-b)*q) : b in [1..(i-1)]];
        newCpq := newCpqs[1];
        b := 1;
        // Choose the cheapest strategy.
        for j in [2..(i-1)] do
            tmpCpq := newCpqs[j];
            if newCpq gt tmpCpq then
                // strict inequality in the condition prefers larger number of isogenies
                newCpq := tmpCpq;
                b := j;
            end if;
        end for;
  	// chosen strategy (m-leave sub-tree on the left, (n-m)-subtree on the right)
       Append(~S,[b] cat S[i-b] cat S[b]);
       Append(~C, newCpqs[b]);
    end for;
    return S[n];
end function;

// Compute large isogeny
// [Adapted from vOW4SIKE]
function PowerOfPrimeIsogeny(a,LEN,C0,STRAT,PTS : inv:=true)
    a24,xp,xq,xpq := Explode(C0);
    X,Z := LADDER_3_pt(a,xp,xq,xpq,a24);
    LENSTRAT := #STRAT;

    POINTS := []; for i:=1 to #PTS do Append(~POINTS,[PTS[i],1]); end for; // Put x into X,Z

    pts := [];
    index := 0;
    ii := 1;

    // Treat the first row separatelyfor i:=1 to #P do Append(~PTS,[P[i],1]); end for;
    // multiply (X:Z) until it has order 4, and store intermediate points
    while index lt LENSTRAT do
        Append(~pts, [X,Z,index]);
        m := STRAT[ii]; ii +:= 1;
        X,Z := xDBLe_affineA24(X,Z,a24,m); // I changed 2m to m (expect it to be related to 2e/2)
        index +:= m;
    end while;

    // compute the 4-isogeny based on kernel (X:Z)
    // evaluate the 4-isogeny at every point in pts
    // map along the points in POINTS
    if X eq Z then
        for i:=1 to #pts do
            pts[i][1],pts[i][2] := EvalFourIsogenyWithKernelXeqZ(a24,pts[i][1],pts[i][2]);
        end for;
	for i:=1 to #POINTS do
	    POINTS[i][1],POINTS[i][2] := EvalFourIsogenyWithKernelXeqZ(a24,POINTS[i][1],POINTS[i][2]);
	end for;
        A24,C24 := GetFourIsogenyWithKernelXeqZ(a24);
    elif X eq -Z then
        for i:=1 to #pts do
            pts[i][1],pts[i][2] := EvalFourIsogenyWithKernelXeqMinZ(a24,pts[i][1],pts[i][2]);
        end for;
	for i:=1 to #POINTS do
            POINTS[i][1],POINTS[i][2] := EvalFourIsogenyWithKernelXeqMinZ(a24,POINTS[i][1],POINTS[i][2]);
        end for;
        A24,C24 := GetFourIsogenyWithKernelXeqMinZ(a24);
    else
        A24,C24,K1,K2,K3 := GetFourIsogenyWithKernelXneZ(X,Z);
        for i:=1 to #pts do
            pts[i][1],pts[i][2] := EvalFourIsogenyWithKernelXneZ(pts[i][1],pts[i][2],K1,K2,K3);
        end for;
	for i:=1 to #POINTS do
            POINTS[i][1],POINTS[i][2] := EvalFourIsogenyWithKernelXneZ(POINTS[i][1],POINTS[i][2],K1,K2,K3);
        end for;
    end if;

    if LENSTRAT gt 0 then // For very small e

        // R becomes the last point in pts and then pts is pruned
        X := pts[#pts][1];
        Z := pts[#pts][2];
        index := Integers()!pts[#pts][3];

        Prune(~pts); //remove last point

        // Alice's main loop
        for row := 2 to LENSTRAT do
            // multiply (X:Z) until it has order 4, and store intermediate points
            while index lt (LENSTRAT + 1 - row) do
                Append(~pts, [X,Z,index]);
                m := STRAT[ii]; ii +:= 1;
                X,Z := xDBLe(X,Z,A24,C24,2*m);
                index +:= m;
            end while;

            // compute the 4-isogeny based on kernel (X:Z)
            A24,C24,K1,K2,K3 := GetFourIsogenyWithKernelXneZ(X,Z);

            // evaluate the 4-isogeny at every point in pts
            for i:=1 to #pts do
                pts[i][1],pts[i][2] := EvalFourIsogenyWithKernelXneZ(pts[i][1],pts[i][2],K1,K2,K3);
            end for;
	    for i:=1 to #POINTS do
                POINTS[i][1],POINTS[i][2] := EvalFourIsogenyWithKernelXneZ(POINTS[i][1],POINTS[i][2],K1,K2,K3);
            end for;


            // R becomes the last point in pts and then pts is pruned
            X := pts[#pts][1];
            Z := pts[#pts][2];
            index := Integers()!pts[#pts][3];

            Prune(~pts);
        end for;

        // compute the last 4-isogeny and map along the points
	A24,C24,K1,K2,K3 := GetFourIsogenyWithKernelXneZ(X,Z);
	for i:=1 to #POINTS do
            POINTS[i][1],POINTS[i][2] := EvalFourIsogenyWithKernelXneZ(POINTS[i][1],POINTS[i][2],K1,K2,K3);
        end for;

    end if;

    A24 := A24 + A24;
    A24 := A24 - C24;
    A24 := A24 + A24;

    //j := j_inv(A24,C24);

    if inv then
        a24 := A24/C24; //Extra inversion needed in delegation, skip this for local runtimes
	Cset := [a24];
    else
	Cset := [A24,C24];
    end if;

    for i:=1 to #POINTS do
	Append(~Cset,POINTS[i][1]/POINTS[i][2]);
    end for;

    return Cset;
end function;


// Initialize parameters
function InitParams(a0,e2,e3,f)
    p := 2^e2*3^e3*f - 1;
    _<x> := PolynomialRing(GF(p));
    Fp2<i> := ext< GF(p) | x^2 + 1>;
    a0 := Fp2!a0;

    a24:= (a0+2)/4;
    E0 := EllipticCurve([Fp2|0,a0,0,1,0]);

    //P2
    repeat P20 := f*Random(E0); P2 := 3^e3*P20; P22 := 2^(e2-1)*P2;
    until P22 ne E0!0 and P22 ne E0![0,0] and 2*P22 eq E0!0; //P1 full order and not (0:0:1)
    //Q2
    repeat Q20 := f*Random(E0); Q2 := 3^e3*Q20; Q22 := 2^(e2-1)*Q2;
    until Q22 eq E0![0,0]; //Q1 not (0:0:1)
    //P3
    repeat P30 := f*Random(E0); P3 := 2^e2*P30; P32 := 3^(e3-1)*P3;
    until P32 ne E0!0; // and P32 ne E0![0,0] and 2*P32 eq E0!0; //P1 full order and not (0:0:1)
    //Q3
    repeat Q30 := f*Random(E0); Q3 := 2^e2*Q30; Q32 := 3^(e3-1)*Q3;
    until Q32 ne E0!0; //Q1 not (0:0:1)

    xP2 := P2[1]; xQ2 := Q2[1]; xPQ2 := (P2-Q2)[1];
    xP3 := P3[1]; xQ3 := Q3[1]; xPQ3 := (P3-Q3)[1];

    return Fp2,[a24,xP2,xQ2,xPQ2], [a24,xP3,xQ3,xPQ3];
end function;









// ----------------------------------- //
// ---- Cryptographic subroutines ---- //
// ----------------------------------- //

function Pseudo_B_Isogeny(b,e2,e3,f,C2,strat2,Points)
    // 3-isogenies are not implemented in vOW4SIKE; for benchmarking purposes,
    //  we simulate them using 2-isogenies. A public-key is simulated by simply
    //  generating a random elliptic curve (b being random)
    CB := PowerOfPrimeIsogeny(b,e2,C2,strat2,Points);
    a24B := CB[1]; Points := CB[2..#CB];
    EB := EllipticCurveFromjInvariant(j_inv(a24B,1));
    // generate 2^e2-torsion basis
    repeat P2 := 3^e3*f*Random(EB); until 2^(e2-1)*P2 ne EB!0;
    repeat Q2 := 3^e3*f*Random(EB); until 2^(e2-1)*Q2 ne EB!0;
    // return results
    return [a24B,P2[1],Q2[1],(P2-Q2)[1]],Points;
end function;

function PublicKey_Alice(a,e2,C2,C3,strat2 : inv:=true)
    _,xP3,xQ3,xPQ3 := Explode(C3);
    return PowerOfPrimeIsogeny(a,e2,C2,strat2,[xP3,xQ3,xPQ3] : inv:=inv);
end function;

function PseudoPublicKey_Bob(b,e2,e3,f,C2,strat2)
    CB,_ := Pseudo_B_Isogeny(b,e2,e3,f,C2,strat2,[]);
    return CB;
end function;

function SharedCurve(a,e2,CB,strat2 : inv:=true) //for Alice
    CB_A := PowerOfPrimeIsogeny(a,e2,CB,strat2,[] : inv:=inv);
    return CB_A;
end function;

function ZKPI(xPA,e2,C2,CA,strat2 : inv:=true)
    b := Integers()!Random(Integers(2^e2));
    CB:= PowerOfPrimeIsogeny(b,e2,C2,strat2,[xPA] : inv:=inv);
    CAB:=PowerOfPrimeIsogeny(b,e2,CA,strat2,[]);
    return CB,CAB;
end function;






// ----------------------------------- //
// ---- Timing local computations ---- //
// ----------------------------------- //




function time_PublicKey_local(a0,e2,e3,f)
    strat2 := OptStrat(e2,costS,costI);
    Fp2<i>,C2,C3 := InitParams(a0,e2,e3,f);

    a := Integers()!Random(Integers(2^e2));
    t := ClockCycles();
        _ := PublicKey_Alice(a,e2,C2,C3,strat2 : inv:=false);
    T := ClockCycles()-t;
    return T;
end function;

function time_SharedCurve_local(a0,e2,e3,f)
    strat2 := OptStrat(e2,costS,costI);
    Fp2<i>,C2,C3 := InitParams(a0,e2,e3,f);

    a := Integers()!Random(Integers(2^e2));
    b := Integers()!Random(Integers(2^e2));
    CB := PseudoPublicKey_Bob(b,e2,e3,f,C2,strat2);
    t := ClockCycles();
	_ := SharedCurve(a,e2,CB,strat2 : inv:=false);
    T := ClockCycles()-t;
    return T;
end function;

function time_ZKPI_local(a0,e2,e3,f)
    strat2 := OptStrat(e2,costS,costI);
    Fp2<i>,C2,C3 := InitParams(a0,e2,e3,f);
    E0 := EllipticCurve([Fp2|0,a0,0,1,0]);

    repeat A := 2^e2*f*Random(E0); until 3^(e3-1)*A ne E0!0; xPA := A[1];
    a := Integers()!Random(Integers(2^e2));
    CA := PseudoPublicKey_Bob(a,e2,e3,f,C2,strat2);
    t := ClockCycles();
        _,_ := ZKPI(xPA,e2,C2,CA,strat2 : inv:=false);
    T := ClockCycles()-t;
    return T;
end function;


