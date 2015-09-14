function E0 = Ekep(M,e)
% eccentric anomaly (by newton-raphson)
% M MUST BE IN RADIANS!
tol = 1e-9; % tolerance
E0 = M + e * sin(M);
de = 1;
while abs(de) > tol
    dm = M - (E0 - e * sin(E0));
    de = dm / (1 - e * cos(E0));
    E1 = E0 + de;
    E0 = E1;
end
return