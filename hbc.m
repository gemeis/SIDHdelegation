// ---- ISO - HBC ---- //

function ISO_HBC_Pre(C2)
    Points := C2[2..#C2];
    return Points;
end function;

function ISO_HBC_U(a,e2,C2,strat2,POINTS)
    CA := PowerOfPrimeIsogeny(a,e2,C2,strat2,POINTS);
    return CA;
end function;

function ISO_HBC_Post(secret,CA)
    a24A,xPA,xQA,xPQA := Explode(CA);
    result := LADDER_3_pt(secret,xPA,xQA,xPQA,a24A);
    return result;
end function;

function time_ISO_HBC(a0,e2,e3,f)
    strat2 := OptStrat(e2 div 2,costS,costI);
    _,C2,_ := InitParams(a0,e2,e3,f);

    a := Integers()!Random(Integers(2^e2));
    s := Integers()!Random(Integers(2^e2));

            t := ClockCycles();
    POINTS := ISO_HBC_Pre(C2);
            T := ClockCycles()-t;
    CA:= ISO_HBC_U(a,e2,C2,strat2,POINTS);
            t := ClockCycles();
    _ := ISO_HBC_Post(s,CA);
            T +:= ClockCycles()-t;

    return T;
end function;









// ---- ISODET HBC ---- //

function ISODET_HBC_Pre(e2,C2,C3)
    _,xP2,xQ2,xPQ2 := Explode(C2);
    _,xP3,xQ3,xPQ3 := Explode(C3);
    Points := [xP2,xQ2,xPQ2,xP3,xQ3,xPQ3];
    k := Integers()!Random(Integers(2^e2));
    return k,Points;
end function;

function ISODET_HBC_U1(k,e2,e3,f,C2,strat2,POINTS)
    CK,POINTS_B := Pseudo_B_Isogeny(k,e2,e3,f,C2,strat2,POINTS[4..#POINTS]);
    return CK,POINTS_B;
end function;

function ISODET_HBC_12(a,CK)
    aK,xP2,xQ2,xPQ2 := Explode(CK);
    A := LADDER_3_pt(a,xP2,xQ2,xPQ2,aK);
    return A;
end function;

function ISODET_HBC_U2(a,e2,CK,strat2,POINTS)
    CKA := PowerOfPrimeIsogeny(a,e2,CK,strat2,POINTS); //here should be A directly instead of a
    return CKA;
end function;

function ISODET_HBC_23(k,CKA)
    a24KA,xPK,xQK,xPQK := Explode(CKA);
    K := LADDER_3_pt(k,xPK,xQK,xPQK,a24KA);
    return K;
end function;

function ISODET_HBC_U3(k,e2,CKA,strat2)
    CA := PowerOfPrimeIsogeny(k,e2,CKA,strat2,[]);
    return CA;
end function;









// ---- HBC full Routines ---- //

function time_ISODET_HBC(a0,e2,e3,f)
    strat2 := OptStrat(e2 div 2,costS,costI);
    Fp2<i>,C2,C3 := InitParams(a0,e2,e3,f);

    a := Integers()!Random(Integers(2^e2));

                t := ClockCycles();
    k,POINTS := ISODET_HBC_Pre(e2,C2,C3);
                T := ClockCycles()-t;
    CK,POINTS := ISODET_HBC_U1(k,e2,e3,f,C2,strat2,POINTS);
                t := ClockCycles();
    _ := ISODET_HBC_12(a,CK);
                T +:= ClockCycles()-t;
    CKA := ISODET_HBC_U2(a,e2,CK,strat2,POINTS);
                t := ClockCycles();
    K := ISODET_HBC_23(k,CKA);
                T +:= ClockCycles()-t;
    CA := ISODET_HBC_U3(k,e2,CKA,strat2);

    return T;
end function;

function time_PUBKEY_HBC(a0,e2,e3,f) //same as IsoDet except that I need to map 2 more points!
    strat2 := OptStrat(e2 div 2,costS,costI);
    Fp2<i>,C2,C3 := InitParams(a0,e2,e3,f);

    a := Integers()!Random(Integers(2^e2));

    t := ClockCycles();
        k,POINTS := ISODET_HBC_Pre(e2,C2,C3);
        for i:=2 to #C2 do Append(~POINTS,C2[i]); end for;
    T := ClockCycles()-t;
        CK,POINTS := ISODET_HBC_U1(k,e2,e3,f,C2,strat2,POINTS);
    t := ClockCycles();
        _ := ISODET_HBC_12(a,CK);
    T +:= ClockCycles()-t;
        CKA := ISODET_HBC_U2(a,e2,CK,strat2,POINTS);
    t := ClockCycles();
        K := ISODET_HBC_23(k,CKA);
    T +:= ClockCycles()-t;
        CA := ISODET_HBC_U3(k,e2,CKA,strat2);
    return T;
end function;

function time_ZKPI_HBC(a0,e2,e3,f)
    strat2 := OptStrat(e2 div 2,costS,costI);
    _,C2,_ := InitParams(a0,e2,e3,f);

    a := Integers()!Random(Integers(2^e2));
    CA := PseudoPublicKey_Bob(a,e2,e3,f,C2,strat2);

    t := ClockCycles();
        b := Integers()!Random(Integers(2^e2));
        Points := ISO_HBC_Pre(C2);
    T := ClockCycles()-t;
        CB := ISO_HBC_U(b,e2,C2,strat2,Points);
        CAB:= ISO_HBC_U(b,e2,CA,strat2,[]);
    t := ClockCycles();
        _ := ISO_HBC_Post(a,CB);
    T +:= ClockCycles()-t;

    return T;
end function;

function time_SIDH_HBC(a0,e2,e3,f)
    strat2 := OptStrat(e2 div 2,costS,costI);
    _,C2,_ := InitParams(a0,e2,e3,f);
    b := Integers()!Random(Integers(2^e2));
    CB := PseudoPublicKey_Bob(b,e2,e3,f,C2,strat2);

    T := time_PUBKEY_HBC(a0,e2,e3,f);
    t := ClockCycles();
        a := Integers()!Random(Integers(2^e2));
        EAB := PowerOfPrimeIsogeny(a,e2,CB,strat2,[]);
    T +:= ClockCycles()-t;
    return T;
end function;

