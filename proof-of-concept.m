function gensmall(t,n) // generate small scalars
    res := [Integers()!Random(Integers(2^(t-1))) : i in [1..n]];
    return res;
end function;

function generators(E,n) // create basis of order n on E
    N := #GroundField(BaseField(E))+1;
    error if not IsDivisibleBy(N,n), "no torsion group of given order";
    repeat P := (N div n)*Random(E); until Order(P) eq n;
    repeat Q := (N div n)*Random(E);
    until (Order(Q) eq n) and (Order(WeilPairing(P,Q,n)) eq n);
    return P,Q;
end function;

function isog(E,a,P,Q,L) //compute isogeny E/<P+aQ> and push through points in L
    Fq<x> := PolynomialRing(BaseField(E));
    K := P+a*Q;
    Fact := Factorisation(Order(K));
    assert #Fact eq 1;
    l := Fact[1][1]; d := Fact[1][2];
    for j:=1 to d do
        Ker := l^(d-j)*K; // compute kernel generator
        //assert Order(K) eq l; // comment to speed up
        E,phi := IsogenyFromKernel(E, &*{(x-(n*Ker)[1]) : n in [1..l-1]});
        L := [phi(Q) : Q in L];
        K := phi(K);
    end for;
    return E,L;
end function;

// verifying routines
function verify_Iso_HBC(E,k,P3,Q3,a,P2,Q2) //kernel P3+k*Q3, push through A=P2+a*Q2
    A := P2+a*Q2;
    L := [P2,Q2,A];
    EK,LK := isog(E,k,P3,Q3,L);
    P2K := LK[1]; Q2K := LK[2]; AK_local := LK[3];
    return AK_local eq P2K+a*Q2K;
end function;

function verify_Iso_OMTUP(E,k,P3,Q3,a,P2,Q2,e2) //kernel P3+k*Q3, push through A=P2+a*Q2
    N := 4;
    t := 3;
    c := gensmall(t,N); //generate N small integers
    d := gensmall(t,N); //idem
    r0:= Random(Integers(2^e2));
    s0:= Random(Integers(2^e2));
    s := [Random(Integers(2^e2)) : i in [1..N-1]];
    r := [-s[i]+c[i]*s0+d[i]*r0 : i in [1..N-1]];
    sigma := 0; for i:=1 to N-1 do sigma +:= s[i]+r[i]; end for;
    gamma := Integers(2^e2)!3; //for l=2; otherwise change appropriately
    Append(~s,1/gamma*(d[N]*r0+c[N]*s0+sigma-Integers()!a)); //s[N]
    Append(~r,-s[N]+c[N]*s0+d[N]*r0); //r[N]

    A := P2+a*Q2;
    L := [P2,Q2,A];
    EK,LK := isog(E,k,P3,Q3,L);
    P2K := LK[1]; Q2K := LK[2]; AK_local := LK[3];

    AK := Integers()!r[N]*Q2K-(Integers()!gamma-1)*Integers()!s[N]*Q2K; //compute map of A
    for i:=1 to N-1 do AK+:=Integers()!s[i]*Q2K+Integers()!r[i]*Q2K; end for;
    AK +:= P2K;

    return AK eq AK_local;
end function;

function verify_IsoDetour(E,a,P2,Q2,k,P3,Q3,X,e2,e3,e5) // pushing through X via tau3,tau2,tau3
    // first isogeny (kappa)
    L0 := [P2,Q2,X,P3,Q3];
    EK,LK := isog(E,k,P3,Q3,L0);
    P2K := LK[1]; Q2K := LK[2]; XK := LK[3]; P3K := LK[4]; Q3K := LK[5];

    // generate "canonical basis" on EK
    //repeat // some SK might be ill-formed
    SK := 2^e2*5^e5*Random(EK);
    RK := P3K-k*SK;
    //until Order(RK+k*SK) eq Order(P3+k*Q3);

    // second isogeny (alpha')
    L1 := [RK,SK,XK];
    EAK,LAK := isog(EK,a,P2K,Q2K,L1);
    RAK := LAK[1]; SAK := LAK[2]; XAK := LAK[3];

    // third isogeny (kappa^')
    L2 := [InverseMod(3^e3,5^e5)*XAK];
    EA,LA := isog(EAK,k,RAK,SAK,L2);
    XA := LA[1];

    // local version
    EA_local,LA_local := isog(E,a,P2,Q2,[X]);
    XA_local := LA_local[1];

    iso := Isomorphism(EA,EA_local);
    return iso(XA) eq XA_local;
end function;

procedure run()
    // Testing Iso
    e2 := 216; e3 := 137;
    p := 2^e2*3^e3-1;
    E := EllipticCurve([GF(p,2)!1,0]);

    P2,Q2 := generators(E,2^e2);
    repeat a := Integers()!Random(Integers(2^e2)); until (Gcd(a,2) eq 1);
    P3,Q3 := generators(E,3^e3);
    repeat k := Integers()!Random(Integers(3^e3)); until (Gcd(k,3) eq 1);

    print "";
    print "Iso_HBC",verify_Iso_HBC(E,k,P3,Q3,a,P2,Q2);
    print "Iso_OMTUP",verify_Iso_OMTUP(E,k,P3,Q3,a,P2,Q2,e2);

    // Testing IsoDet
    //e2 := 207; e3 := 131; e5 := 93; //slow
    e2 := 63; e3 := 41; e5 := 27; //faster for testing
    p := 2^e2*3^e3*5^e5-1;
    E := EllipticCurve([GF(p,2)!1,0]);

    X := 2^e2*3^e3*Random(E); // random 5^e5-torsion point

    P2,Q2 := generators(E,2^e2);
    repeat a := Integers()!Random(Integers(2^e2)); until (Gcd(a,2) eq 1);
    P3,Q3 := generators(E,3^e3);
    repeat k := Integers()!Random(Integers(3^e3)); until (Gcd(k,3) eq 1);

    print "IsoDetour", verify_IsoDetour(E,a,P2,Q2,k,P3,Q3,X,e2,e3,e5);
end procedure;

run();
