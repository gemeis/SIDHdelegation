function EADD(P1,P2,a,d) //9m
    X1,Y1,T1,Z1 := Explode(P1);
    X2,Y2,T2,Z2 := Explode(P2);
    A := X1*X2;
    B := Y1*Y2;
    C := T1*d*T2;
    D := Z1*Z2;
    E := (X1+Y1)*(X2+Y2)-A-B;
    F := D-C;
    G := D+C;
    H := B-a*A;
    X3 := E*F;
    Y3 := G*H;
    T3 := E*H;
    Z3 := F*G;
    return [X3,Y3,T3,Z3];
end function;

function EDBL(P,a,d)
    X,Y,T,Z := Explode(P);
    A := X^2;
    B := Y^2;
    C := 2*Z^2;
    D := a*A;
    E := (X+Y)^2-A-B;
    G := D+B;
    F := G-C;
    H := D-B;
    X3 := E*F;
    Y3 := G*H;
    T3 := E*H;
    Z3 := F*G;
    return [X3,Y3,T3,Z3];
end function;

function EDBLe(m,P,a,d)
    bits:=IntegerToSequence(m,2);
    Q := [0,0,0,0];
    for i:=1 to #bits do
        if bits[i] eq 0 then P:=EDBL(P,a,d); end if;
        if bits[i] eq 1 then
	    if Q[4] eq 0 then Q:=P;
	    else Q := EADD(Q,P,a,d);
	    end if;
	end if;
    end for;
    return Q;
end function;

function Mont2Edwards(a24,xPoints) // translate montgomery to Edwards
    a0 := 4*a24-2;
    c := a0+2; d := a0-2;
    E := EllipticCurve([0,a0,0,1,0]);
    PointsEd := [];
    for x in xPoints do
        bool,P := IsPoint(E,x);
        if bool then
            xP:=P[1]; yP:=P[2];
        else // randomize for benchmarking purposes
            P := Random(E); xP := P[1]; yP := P[2];
        end if;
        xPE := xP*(xP+1);
        yPE := yP*(xP-1);
        Append(~PointsEd,[xPE,yPE,xPE*yPE,1]);
    end for;
    return c,d,PointsEd;
end function;
