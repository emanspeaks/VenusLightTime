close all
clear
clc

% utility functions
deg2rad = @(x)(x*pi/180);
rad2deg = @(x)(x*180/pi);
greg2jd = @(m,d,y,ut)(367*y-floor((7*(y+floor((m+9)/12)))/4)+floor((275*m)/9)+d+1721013.5+ut/24-0.5*sign(100*y+m-190002.5)+0.5);
pm180 = @(x)(mod(x+180,360)-180);
pmpi = @(x)(mod(x+pi,2*pi)-pi);

% constants
au = 1.49597870691e8; % astronomical unit (km)
mu = 1.32712440018e11; % heliocentric gravity parameter (km^3/s^2)

% planetary data for 1800 AD - 2050 AD
% from http://ssd.jpl.nasa.gov/txt/p_elem_t1.txt
data = [0.38709927	,	0.20563593	,	7.00497902	,	252.2503235	,	77.45779628	,	48.33076593	;
0.00000037	,	0.00001906	,	-0.00594749	,	149472.6741	,	0.16047689	,	-0.12534081	;
0.72333566	,	0.00677672	,	3.39467605	,	181.9790995	,	131.6024672	,	76.67984255	;
0.0000039	,	-0.00004107	,	-0.0007889	,	58517.81539	,	0.00268329	,	-0.27769418	;
1.00000261	,	0.01671123	,	-0.00001531	,	100.4645717	,	102.9376819	,	0	;
0.00000562	,	-0.00004392	,	-0.01294668	,	35999.37245	,	0.32327364	,	0	;
1.52371034	,	0.0933941	,	1.84969142	,	-4.55343205	,	-23.94362959	,	49.55953891	;
0.00001847	,	0.00007882	,	-0.00813131	,	19140.30268	,	0.44441088	,	-0.29257343	;
5.202887	,	0.04838624	,	1.30439695	,	34.39644051	,	14.72847983	,	100.4739091	;
-0.00011607	,	-0.00013253	,	-0.00183714	,	3034.746128	,	0.21252668	,	0.20469106	;
9.53667594	,	0.05386179	,	2.48599187	,	49.95424423	,	92.59887831	,	113.6624245	;
-0.0012506	,	-0.00050991	,	0.00193609	,	1222.493622	,	-0.41897216	,	-0.28867794	;
19.18916464	,	0.04725744	,	0.77263783	,	313.2381045	,	170.9542763	,	74.01692503	;
-0.00196176	,	-0.00004397	,	-0.00242939	,	428.4820279	,	0.40805281	,	0.04240589	;
30.06992276	,	0.00859048	,	1.77004347	,	-55.12002969	,	44.96476227	,	131.7842257	;
0.00026291	,	0.00005105	,	0.00035372	,	218.4594533	,	-0.32241464	,	-0.00508664	;
39.48211675	,	0.2488273	,	17.14001206	,	238.9290383	,	224.0689163	,	110.3039368	;
-0.00031596	,	0.0000517	,	0.00004818	,	145.2078052	,	-0.04062942	,	-0.01183482	];

% functions to parse planetary data, p is planet number (1=mercury, 9=pluto)
tpl = @(jd)((jd - 2451545) / 36525);
plval = @(p,jd,el)(data(1+(p-1)*2,el)+data(p*2,el)*tpl(jd));

apl = @(p,jd)(plval(p,jd,1)); % sma, au
epl = @(p,jd)(plval(p,jd,2)); % ecc
ipl = @(p,jd)(deg2rad(plval(p,jd,3))); % inc, rad
lpl = @(p,jd)(deg2rad(plval(p,jd,4))); % mean longitude, rad
wbpl = @(p,jd)(deg2rad(plval(p,jd,5))); % longitude of peri, rad
opl = @(p,jd)(deg2rad(plval(p,jd,6))); % LAN, rad

wpl = @(p,jd)(wbpl(p,jd) - opl(p,jd)); % argument of peri, rad
Mpl = @(p,jd)(pmpi(lpl(p,jd) - wbpl(p,jd))); % mean anomaly [-pi,pi], rad

% generic orbital mechanics relations
E2nu = @(E,e)(pmpi(2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2)))); % true anomaly from ecc anomaly [-pi,pi], rad
r = @(a,e,nu)(a * (1 - e^2) / (1 + e * cos(nu))); % range, units of a

lvlh2xyz = @(co,so,ci,si,cw,sw)([co*cw-so*ci*sw,-co*sw-so*ci*cw,so*si;so*cw+co*ci*sw,-so*sw+co*ci*cw,-co*si;si*sw,si*cw,ci]);
		

% r1 = Earth on April 5, 1972
jd1 = greg2jd(4,5,1972,0);

a1 = apl(3,jd1)*au;
e1 = epl(3,jd1);
i1 = ipl(3,jd1);
w1 = wpl(3,jd1);
o1 = opl(3,jd1);
M1 = Mpl(3,jd1);

E1 = Ekep(M1,e1);
nu1 = E2nu(E1,e1);
rmag1 = r(a1,e1,nu1);

rhat1 = [cos(nu1);sin(nu1);0];
co1 = cos(o1);
so1 = sin(o1);
ci1 = cos(i1);
si1 = sin(i1);
cw1 = cos(w1);
sw1 = sin(w1);
r1 = rmag1 * lvlh2xyz(co1,so1,ci1,si1,cw1,sw1)*rhat1;

% r2 = Venus on August 23, 1972
jd2 = greg2jd(8,23,1972,0);

a2 = apl(2,jd2)*au;
e2 = epl(2,jd2);
i2 = ipl(2,jd2);
w2 = wpl(2,jd2);
o2 = opl(2,jd2);
M2 = Mpl(2,jd2);

E2 = Ekep(M2,e2);
nu2 = E2nu(E2,e2);
rmag2 = r(a2,e2,nu2);

rhat2 = [cos(nu2);sin(nu2);0];
co2 = cos(o2);
so2 = sin(o2);
ci2 = cos(i2);
si2 = sin(i2);
cw2 = cos(w2);
sw2 = sin(w2);
r2 = rmag2 * lvlh2xyz(co2,so2,ci2,si2,cw2,sw2)*rhat2;

% Earth to Venus arc (r1 to r2)
typenum1 = 1; % 1 = short arc, 2 = long arc
TOF1 = (jd2-jd1)*86400; %s

% define Lambert arc plane
typeid1 = typenum1-1; % dot(r1xr2hat,zhat) < 0 (assumes sign(dot(hhat,zhat)) = 1).  If you specify as typenum, dot product not needed.
type1p1 = (1-2*typeid1);
r1xr2 = cross(r1,r2);
r1xr2hat = r1xr2/norm(r1xr2);
hhat1 = type1p1*r1xr2hat; % set angular momentum direction from type

% compute Lambert arc
dnu1 = atan2(norm(r1xr2),dot(r1,r2))*180/pi; %deg
TA1 = 360*typeid1 + dnu1*type1p1; %deg
[aL1,pL1,typeL1,dE1] = lambert(mu,rmag1,rmag2,TA1,TOF1); %km, km, (string), rad

aL1
pL1
typeL1

n1 = sqrt(mu/aL1^3);
eL1 = sqrt(1-pL1/aL1);
hmag1 = sqrt(mu*pL1); %au^2/s
h1 = hmag1*hhat1;

f1 = 1-aL1/rmag1*(1-cos(dE1));
g1 = TOF1-(dE1-sin(dE1))/n1;
v1 = (r2-f1*r1)/g1;
vmag1 = norm(v1);
gamma1 = sign(dot(r1,v1))*acos(hmag1/rmag1/vmag1);
EL1 = sign(gamma1)*acos((1-rmag1/aL1)/eL1);
ML1 = EL1 - eL1*sin(EL1);
nuL1 = E2nu(EL1,eL1);

t = jd1;
Mx = ML1;
Ex = EL1;
rx = r1;
vx = v1;
while t < jd2
	% step forward one day
    t = t + 1;
    My = Mx + n1*86400;
    [m,d,y,u] = jd2greg(t);
    fprintf('%d/%d/%d %dh UT\n',m,d,y,u)
    
    % solve for new state using f and g relations
    Ey = Ekep(My,eL1);
    dEL = mod(Ey-Ex,2*pi);
    fL = 1-aL1/norm(rx)*(1-cos(dEL));
    gL = 86400-(dEL-sin(dEL))/n1;
    ry = rx*fL+vx*gL;
    fdotL = -sqrt(mu*aL1)/norm(rx)/norm(ry)*sin(dEL);
    gdotL = 1-aL1/norm(ry)*(1-cos(dEL));
    vy = rx*fdotL+vx*gdotL;

    % recompute Earth state
    aE = apl(3,t)*au;
    eE = epl(3,t);
    iE = ipl(3,t);
    wE = wpl(3,t);
    oE = opl(3,t);
    ME = Mpl(3,t);

    EE = Ekep(ME,eE);
    nuE = E2nu(EE,eE);
    rmagE = r(aE,eE,nuE);

    rhatE = [cos(nuE);sin(nuE);0];
    coE = cos(oE);
    soE = sin(oE);
    ciE = cos(iE);
    siE = sin(iE);
    cwE = cos(wE);
    swE = sin(wE);
    rE = rmagE * lvlh2xyz(coE,soE,ciE,siE,cwE,swE)*rhatE;

    % compute distance to Earth center in km and light time
    rad2deg(nuE)
    rho = ry - rE;
    rhoKm = norm(rho);
    rhoER = rhoKm/6378;
    rhoLight = rhoKm/299792.458;
    plotdata(t-jd1,1) = t;
    plotdata(t-jd1,2) = rhoKm;
    plotdata(t-jd1,3) = rhoER;
    plotdata(t-jd1,4) = rhoLight;

    % update epoch state for loop
    rx = ry;
    vx = vy;
    Mx = My;
    Ex = Ey;
end

nu = [0:0.1:360];
r_e = a1.*(1-e1.^2)./(1+e1.*cosd(nu));
r_v = a2.*(1-e2.^2)./(1+e2.*cosd(nu));
plot(r_e.*cos(deg2rad(nu)-w1-o1),r_e.*sin(deg2rad(nu)-w1-o1))
axis square
hold on
plot(r_v.*cos(deg2rad(nu)-w2-o2),r_v.*sin(deg2rad(nu)-w2-o2))
r_t = aL1.*(1-eL1.^2)./(1+eL1.*cosd(nu));
plot(r_t.*cos(deg2rad(nu)-nuL1+nu1-w1-o1),r_t.*sin(deg2rad(nu)-nuL1+nu1-w1-o1))
hold off

figure
plot(plotdata(:,4))