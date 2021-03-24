load 'iso.m';
load 'hbc.m';
load 'omtup.m';


procedure benchmark(a0,e2,e3,f,n,N)
    Fp2<i>,_,_ := InitParams(a0,e2,e3,f);

    // -------------------------------------------
    printf "Finite field multiplication:           ";
    T := 0; for i:=1 to n do
        a := Random(Fp2); b := Random(Fp2);
        t := ClockCycles(); _ := a*b; T +:= ClockCycles()-t;
    end for;
    m := T div n; print m,1;

    // -------------------------------------------
    printf "Finite field inversion:               ";
    T := 0; for i:=1 to n do
        a := Random(Fp2);
        t := ClockCycles(); _ := 1/a; T +:= ClockCycles()-t;
    end for;
    inv := T div n;
    print inv, inv div m;


    print "";

    // -------------------------------------------
    printf "Public key generation:            ";
    T := 0; for i:=1 to N do T+:= time_PublicKey_local(a0,e2,e3,f); end for;
    print T div N,T div (N*m);

    // -------------------------------------------
    printf "Shared Curve generation:          ";
    T := 0; for i:=1 to N do T+:= time_SharedCurve_local(a0,e2,e3,f); end for;
    print T div N, T div (N*m);

    // -------------------------------------------
    printf "ZKPI-Prover:                      ";
    T := 0; for i:=1 to N do T+:= time_ZKPI_local(a0,e2,e3,f); end for;
    print T div N, T div (N*m);


    print "";

    // -------------------------------------------
    printf "Delegated Iso_HBC:                 ";
    T := 0; for i:=1 to N do T +:= time_ISO_HBC(a0,e2,e3,f); end for;
    print T div N, T div (N*m);

    // -------------------------------------------
    printf "Delegated IsoDet_HBC:              ";
    T := 0; for i:=1 to N do T +:=time_ISODET_HBC(a0,e2,e3,f); end for;
    print T div N, T div (N*m);

    // -------------------------------------------
    printf "Delegated PubKey_HBC:              ";
    T := 0; for i:=1 to N do T +:=time_PUBKEY_HBC(a0,e2,e3,f); end for;
    print T div N, T div (N*m);

    // -------------------------------------------
    printf "Delegated ZKPI_HBC:                ";
    T := 0; for i:=1 to N do T +:= time_ZKPI_HBC(a0,e2,e3,f); end for;
    print T div N, T div (N*m);

    // -------------------------------------------
    printf "Delegated SIDH_HBC:               ";
    T := 0; for i:=1 to N do T +:= time_SIDH_HBC(a0,e2,e3,f); end for;
    print T div N, T div (N*m);



    print "";

    // -------------------------------------------
    printf "Delegated Iso_OMTUP (t=16,N=1):     ";
    T := 0; for i:=1 to N do T +:= time_ISO_OMTUP(a0,e2,e3,f,16,1); end for;
    print T div N, T div (N*m);

    // -------------------------------------------
    printf "Delegated Iso_OMTUP(t=lam/N,N=4):  ";
    T := 0; for i:=1 to N do T +:= time_ISO_OMTUP(a0,e2,e3,f,e2 div 8,4); end for;
    print T div N, T div (N*m);

    // -------------------------------------------
    printf "Delegated IsoDet_OMTUP:            ";
    T := 0; for i:=1 to N do T +:= time_ISODET_OMTUP(a0,e2,e3,f); end for;
    print T div N, T div (N*m);

    // -------------------------------------------
    printf "Delegated PubKey_OMTUP:            ";
    T := 0; for i:=1 to N do T +:= time_PUBKEY_OMTUP(a0,e2,e3,f); end for;
    print T div N, T div (N*m);

    // -------------------------------------------
    printf "Delegated ZKPI_OMTUP:              ";
    T := 0; for i:=1 to N do T +:= time_ZKPI_OMTUP(a0,e2,e3,f); end for;
    print T div N, T div (N*m);

    // -------------------------------------------
    printf "Delegated SIDH_OMTUP:             ";
    T := 0; for i:=1 to N do T +:= time_SIDH_OMTUP(a0,e2,e3,f); end for;
    print T div N, T div (N*m);
end procedure;


procedure run()
    e2:=216; e3:=137; f:=1; a0:=6;
    L := [<216,137>,<250,159>,<305,192>,<372,239>];
    for l in L do
        e2 := l[1]; e3 := l[2]; f := 1;
        printf("\n\n(e2,e3) = ");
	print e2,e3,"\n";
        benchmark(a0,e2,e3,f,1000,5); //reduce iterations for faster algorithm
    end for;
end procedure;


run();
