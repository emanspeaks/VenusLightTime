function [month day year ut] = jd2greg(input)
jd = input;
x = jd - floor(jd);
if (x>=0.5)
    jd = ceil(jd);
end
l = floor(jd + 68569);
n = floor(4*l/146097);
l = l-floor((146097*n+3)/4);
i = floor(4000*(l+1)/1461001);
l = l-floor(1461*i/4) + 31;
j = floor(80*l/2447);
k = l-floor(2447*j/80);
l = floor(j/11);
j = j + 2 - 12*l;
i = 100*(n-49)+i+l;

year = i;
month = j;
day = k;

ut = x + 0.5;
ut = mod(ut,1)*24;
return