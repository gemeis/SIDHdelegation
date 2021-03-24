load 'edwards.m';

function gensmall(t,e2,n) // generate small scalars
    res := [Integers()!Random(Integers(2^(t-1))) : i in [1..n]]; //don't need to care about negative!
    return res;
end function;


// ---- OPTIMAL MULTIPLICATION PROCEDURE ---- //

// Subroutine to add different sets
function Add_sets(Q,SubBoxes,a,d) //L is the location set
    for B in SubBoxes do
        if B[4] ne 0 then
            if Q[4] eq 0 then Q:=B;
            else Q := EADD(Q,B,a,d);
            end if; end if;
    end for;
    return Q;
end function;


// Subroutine to do many multiplications of the same Q optimally
function Optimal_multiplication(c,Q,a,d,N)

    // Input set c and point Q=[X:Y:T:Z] then returns [c[i]*Q] for i=[1..N]
    // other parameters are the curve a,d and N

    if N eq 1 then return EDBLe(c[1],Q,a,d); end if;

    cbits := [IntegerToSequence(c[i],2) : i in [1..N]]; // put c into bit-sequences
    max := Max([#i : i in cbits]);                      // longest array in cbits
    for i:=1 to N do while #cbits[i] lt max do Append(~cbits[i],0); end while; end for; //add zeros

    dbl := [Q];                            // first element is Q
    for i:=2 to max do                     // Add multiples of Q up to max
        Q2 := EDBL(Q,a,d);                 // so that dbl[i] = 2^(i-1)*Q
        Append(~dbl,Q);
    end for;

    /* We define 2^N-1 "boxes" in the order
        A,B,AB,C,AC,BC,ABC,D,AD,BD,ABD,CD,ACD,BCD,ABCD,E,....
        so that we can map to them using binary strings,
        e.g. A=[1],B=[0,1],AB=[1,1] etc --> ex: ACD=[1,0,1,1], BC=[0,1,1,0] etc*/

    Fp2 := Parent(Q[1]); //needs to be defined so assignment to boxes doesn't throw error
    Boxes := [[Fp2!0,Fp2!0,Fp2!0,Fp2!0] : i in [1..2^N-1]];     // Initialize empty boxes
    for i:=1 to max do                                  	// Add 2^(i-1)*Q to appropriate box
        loc := SequenceToInteger([cbits[j][i] : j in [1..N]],2);//maps e.g. [1,0,1,1] to 13 (ACD)
        if loc ne 0 then
            if Boxes[loc] eq [0,0,0,0] then                 // If empty, initialize
                Boxes[loc][1]:=dbl[i][1]; Boxes[loc][2]:=dbl[i][2];
                Boxes[loc][3]:=dbl[i][3]; Boxes[loc][4]:=dbl[i][4];
            else                                        // else Boxes[loc]+:=dbl[i]
                Boxes[loc][1],Boxes[loc][2],Boxes[loc][3],Boxes[loc][4]
                        := Explode(EADD(Boxes[loc],dbl[i],a,d));
            end if;
        end if;
    end for;

    // We define 15 "boxes" in the order
    // A,B,AB,C,AC,BC,ABC,D,AD,BD,ABD,CD,ACD,BCD,ABCD
    // from here on it's specifically for N=4, other cases not implemented yet
    assert N eq 4;
    KA := Add_sets([0,0,0,0],[Boxes[i] : i in [15,7,11,3]],a,d);   // KA = ABCD+ABC+ABD+AB
    KD := Add_sets([0,0,0,0],[Boxes[i] : i in [15,14,13,12]],a,d); // KD = ABCD+BCD+ACD+CD
    A := Add_sets(KA,[Boxes[i] : i in [13,5,9,1]],a,d);  // A  = KA +ACD+AC+AD+A
    B := Add_sets(KA,[Boxes[i] : i in [14,6,10,2]],a,d); // B  = KA +BCD+BC+BD+B
    C := Add_sets(KD,[Boxes[i] : i in [7,6,5,4]],a,d);   // C  = KD +ABC+BC+AC+C
    D := Add_sets(KD,[Boxes[i] : i in [11,10,9,8]],a,d); // D  = KD +DBA+DB+DA+D

    return [A,B,C,D];
end function;






// ---- ISO - OMTUP ---- //

function ISO_OMTUP_Pre(a,e2,t,N)
    c := gensmall(t,e2,N);
    d := gensmall(t,e2,N);
    r0:= Random(Integers(2^e2));
    s0:= Random(Integers(2^e2));
    s := [Random(Integers(2^e2)) : i in [1..N-1]];
    r := [-s[i]+c[i]*s0+d[i]*r0 : i in [1..N-1]];
    sigma := 0; for i:=1 to N-1 do sigma +:= s[i]+r[i]; end for;
    gamma := Integers(2^e2)!3; //for l=2; otherwise change appropriately
    Append(~s,1/gamma*(d[N]*r0+c[N]*s0+sigma-Integers()!a)); //s[N]
    Append(~r,-s[N]+c[N]*s0+d[N]*r0);
    return c,d,r,s,r0,s0;
end function;

function ISO_OMTUP_U(a,e2,N,rs0,rs,C2,strat2,Points)
    CA := PowerOfPrimeIsogeny(a,e2,C2,strat2,C2[2..#C2] cat Points);
    a24 := CA[1]; xP := CA[2]; xQ := CA[3];
    Points := CA[2..#CA];

    X2Q,Z2Q := xDBL(xQ,1,a24,1); x2Q := X2Q/Z2Q;
    Xrs0Q,Zrs0Q := LADDER_3_pt(Integers()!rs0,xQ,x2Q,xQ,a24); //init: Q,2*Q,Q
    xrs0Q := Xrs0Q/Zrs0Q;
    xrsQ := [];
    for i:=1 to N do
        X,Z := LADDER_3_pt(Integers()!rs[i],xQ,x2Q,xQ,a24);
        Append(~xrsQ,X/Z);
    end for;
    // Translate to Edwards
    _,_,rs0QE := Mont2Edwards(a24,[xrs0Q]); rs0QE := rs0QE[1];
    c,d, rsQE := Mont2Edwards(a24,xrsQ);
    return a24,c,d,rs0QE,rsQE,CA[2..4] cat Points;
end function;

function ISO_OMTUP_Post(cE,cE2,dE,dE2,r0Q,s0Q,rQ,sQ,c,d)
    // E1,E2 = [c,d];
    if ((cE ne cE2) or (dE ne dE2)) then error "Elliptic Curves not equal"; end if;
    N := #rQ;
    LHS := []; sigma := [0,0,0,0];
    for i:=1 to N do
        P := EADD(sQ[i],rQ[i],cE,dE); Append(~LHS,P);
        if sigma[4] eq 0 then sigma:=P; else sigma:=EADD(sigma,P,cE,dE); end if;
    end for;
    // change to optimal
    cs0Q := Optimal_multiplication(c,s0Q,cE,dE,N);
    dr0Q := Optimal_multiplication(d,r0Q,cE,dE,N);
    if N eq 1 then RHS := [EADD(cs0Q,dr0Q,cE,dE)];
    else
        RHS := [];
        for i:=1 to N do
            Append(~RHS,EADD(cs0Q[i],dr0Q[i],cE,dE));
        end for;
    end if;
    // if LHS eq RHS then print "EQ equal"; else print "EQ not equal"; end if;
    // gamma-1=2 for l=2, otherwise change appropriately
    // result := sigma-(gamma-1)*sQ[N]+rQ[N];
    Res := EDBL(sQ[N],cE,dE); //gamma-1
    Res := EADD(Res,rQ[N],cE,dE);
    Res := [-Res[1],Res[2],-Res[3],Res[4]]; //-Res
    Res := EADD(Res,sigma,cE,dE);
    return cE,dE,Res;
end function;


function time_ISO_OMTUP(a0,e2,e3,f,t,N)
    strat2 := OptStrat(e2,costS,costI);
    _,C2,_ := InitParams(a0,e2,e3,f);

    a := Integers()!Random(Integers(2^e2)); //secret
    k := Integers()!Random(Integers(2^e2)); //kernel

    tt:= ClockCycles();
    c,d,r,s,r0,s0 := ISO_OMTUP_Pre(a,e2,t,N);
    T := ClockCycles()-tt;

    _,cE1,dE1,r0Q,rQ,_ := ISO_OMTUP_U(k,e2,N,r0,r,C2,strat2,[]);
    _,cE2,dE2,s0Q,sQ,_ := ISO_OMTUP_U(k,e2,N,s0,s,C2,strat2,[]);

    tt:= ClockCycles();
    aQ := ISO_OMTUP_Post(cE1,cE2,dE1,dE2,r0Q,s0Q,rQ,sQ,c,d);
    T+:= ClockCycles()-tt;

    return T;
end function;




// ---- ISODETOUR - OMTUP ---- //


function time_ISODET_OMTUP(a0,e2,e3,f)
    strat2 := OptStrat(e2,costS,costI);
    Fp2<i>,C2,_ := InitParams(a0,e2,e3,f);

    a := Integers()!Random(Integers(2^e2)); //secret
    // these should be actual generators, but are random for benchmarking purposes
    xPI := Random(Fp2); xQI := Random(Fp2);

    // Isogeny #1
    tt:= ClockCycles();
     kf := Integers()!Random(Integers(2^e2)); // detour kernel
     kb := Integers()!Random(Integers(2^e2)); // detour dual
     c,d,r,s,r0,s0 := ISO_OMTUP_Pre(a,e2,e2 div 8,4); // shrouding a, t=e2/8, N=4
    T := ClockCycles()-tt;

    N := 4;
    a24,cE1,dE1,r0Q,rQ,Points := ISO_OMTUP_U(kf,e2,N,r0,r,C2,strat2,[xPI,xQI]);
     _ ,cE2,dE2,s0Q,sQ,Points := ISO_OMTUP_U(kf,e2,N,s0,s,C2,strat2,[xPI,xQI]);
    CK := [a24] cat Points;

    tt:= ClockCycles();
     aQ := ISO_OMTUP_Post(cE1,cE2,dE1,dE2,r0Q,s0Q,rQ,sQ,c,d);
    T+:= ClockCycles()-tt;


    // Isogeny #2
    tt:= ClockCycles();
     c,d,r,s,r0,s0 := ISO_OMTUP_Pre(kb,e2,e2 div 8,4); // shrouding kb, t=e2/8, N=4
    T+:= ClockCycles()-tt;

    a24,cE1,dE1,r0Q,rQ,Points := ISO_OMTUP_U(a,e2,N,r0,r,CK,strat2,[]);
     _ ,cE2,dE2,s0Q,sQ,Points := ISO_OMTUP_U(a,e2,N,s0,s,CK,strat2,[]);
    CAK := [a24] cat Points;

    tt:= ClockCycles();
     kbQ := ISO_OMTUP_Post(cE1,cE2,dE1,dE2,r0Q,s0Q,rQ,sQ,c,d);
    T+:= ClockCycles()-tt;

    a24,cE1,dE1,r0Q,rQ,Points := ISO_OMTUP_U(kb,e2,N,r0,r,CAK,strat2,[]);
     _ ,cE2,dE2,s0Q,sQ,Points := ISO_OMTUP_U(kb,e2,N,s0,s,CAK,strat2,[]);

    return T;
end function;






// ---- OMTUP full Routines ---- //

function time_PUBKEY_OMTUP(a0,e2,e3,f)
    strat2 := OptStrat(e2,costS,costI);
    Fp2<i>,C2,_ := InitParams(a0,e2,e3,f);

    a := Integers()!Random(Integers(2^e2)); //secret
    // these should be actual generators, but are random for benchmarking purposes
    xPI := Random(Fp2); xQI := Random(Fp2);
    xPB := Random(Fp2); xQB := Random(Fp2);

    // Isogeny #1
    tt:= ClockCycles();
     kf := Integers()!Random(Integers(2^e2)); // detour kernel
     kb := Integers()!Random(Integers(2^e2)); // detour dual
     c,d,r,s,r0,s0 := ISO_OMTUP_Pre(a,e2,e2 div 8,4); // shrouding a, t=e2/8, N=4
    T := ClockCycles()-tt;

    N := 4;
    a24,cE1,dE1,r0Q,rQ,Points := ISO_OMTUP_U(kf,e2,N,r0,r,C2,strat2,[xPI,xQI,xPB,xQB]);
     _ ,cE2,dE2,s0Q,sQ,Points := ISO_OMTUP_U(kf,e2,N,s0,s,C2,strat2,[xPI,xQI,xPB,xQB]);
    CK := [a24] cat Points;

    tt:= ClockCycles();
     aQ := ISO_OMTUP_Post(cE1,cE2,dE1,dE2,r0Q,s0Q,rQ,sQ,c,d);
    T+:= ClockCycles()-tt;

    // Isogeny #2
    tt:= ClockCycles();
     c,d,r,s,r0,s0 := ISO_OMTUP_Pre(kb,e2,e2 div 8,4); // shrouding kb, t=e2/8, N=4
    T+:= ClockCycles()-tt;

    a24,cE1,dE1,r0Q,rQ,Points := ISO_OMTUP_U(a,e2,N,r0,r,CK,strat2,[xPB,xQB]);
     _ ,cE2,dE2,s0Q,sQ,Points := ISO_OMTUP_U(a,e2,N,s0,s,CK,strat2,[xPB,xQB]);
    CAK := [a24] cat Points;

    tt:= ClockCycles();
     kbQ := ISO_OMTUP_Post(cE1,cE2,dE1,dE2,r0Q,s0Q,rQ,sQ,c,d);
    T+:= ClockCycles()-tt;


    // Isogeny #3
    a24,cE1,dE1,r0Q,rQ,Points := ISO_OMTUP_U(kb,e2,N,r0,r,CAK,strat2,[xPB,xQB]);
     _ ,cE2,dE2,s0Q,sQ,Points := ISO_OMTUP_U(kb,e2,N,s0,s,CAK,strat2,[xPB,xQB]);

    return T;
end function;




function time_ZKPI_OMTUP(a0,e2,e3,f)
    strat2 := OptStrat(e2,costS,costI);
    _,C2,_ := InitParams(a0,e2,e3,f);

    a := Integers()!Random(Integers(2^e2));
    CA := PseudoPublicKey_Bob(a,e2,e3,f,C2,strat2);

    t := ClockCycles();
        b := Integers()!Random(Integers(2^e2));
	c,d,r,s,r0,s0 := ISO_OMTUP_Pre(a,e2,e2 div 8,4);
    T := ClockCycles()-t;

    N := 4;
    a24,cE1,dE1,r0Q,rQ,Points := ISO_OMTUP_U(b,e2,N,r0,r,C2,strat2,[]);
     _ ,cE2,dE2,s0Q,sQ,Points := ISO_OMTUP_U(b,e2,N,s0,s,C2,strat2,[]);
    _ := PowerOfPrimeIsogeny(b,e2,CA,strat2,[]);
    _ := PowerOfPrimeIsogeny(b,e2,CA,strat2,[]);

    t := ClockCycles();
        aQ := ISO_OMTUP_Post(cE1,cE2,dE1,dE2,r0Q,s0Q,rQ,sQ,c,d);
    T +:= ClockCycles()-t;
    return T;
end function;




function time_SIDH_OMTUP(a0,e2,e3,f)
    strat2 := OptStrat(e2,costS,costI);
    _,C2,_ := InitParams(a0,e2,e3,f);
    b := Integers()!Random(Integers(2^e2));
    CB := PseudoPublicKey_Bob(b,e2,e3,f,C2,strat2);

    T := time_PUBKEY_OMTUP(a0,e2,e3,f);
    t := ClockCycles();
        a := Integers()!Random(Integers(2^e2));
        EAB := PowerOfPrimeIsogeny(a,e2,CB,strat2,[]);
    T +:= ClockCycles()-t;
    return T;
end function;

