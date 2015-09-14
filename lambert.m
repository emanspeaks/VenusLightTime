function [a,p,type,d_anomaly] = lambert(mu,r1,r2,TA,TOF)
% mu and TOF must have same time units
% mu, r1, and r2 must have same distance units

typeid = TA > 180;
typenum = typeid + 1;
type1p = (1-2*typeid);

rho = 360*typeid + TA*type1p; %deg
c = sqrt(r1^2+r2^2-2*r1*r2*cosd(rho));
s = (r1+r2+c)/2;

TOFpara = (s^(3/2)-type1p*(s-c)^(3/2))/3*sqrt(2/mu); %s

switch sign(TOF - TOFpara)
    case -1
        typeconic = 'H';
        
        % alpha equations
        alpha0 = @(a)(2*asinh(sqrt(s/2/abs(a))));
        dalpha0da = @(a)(sqrt(s/(2*abs(a)+s)/a^2));
        alpha = @(a)(alpha0(a));
        dalphada = @(a)(dalpha0da(a));
        
        % beta equations
        beta0 = @(a)(2*asinh(sqrt((s-c)/2/abs(a))));
        dbeta0da = @(a)(sqrt((s-c)/(2*abs(a)+s-c)/a^2));
        beta = @(a)(type1p*beta0(a));
        dbetada = @(a)(type1p*dbeta0da(a));
        
        % mean motion equations
        n = @(a)(sqrt(mu/abs(a)^3));
        dnda = @(a)(sqrt(mu)*3/2*abs(a)^(-5/2));
        
        % Lambert equations
        Q = @(E)(sinh(E)-E);
        dQdE = @(E)(cosh(E)-1);
        dM = @(a)(Q(alpha(a))-Q(beta(a)));
        dMda = @(a)(dQdE(alpha(a))*dalphada(a)-...
            dQdE(beta(a))*dbetada(a));
        K = @(a)(dM(a)-n(a)*TOF);
        Kdot = @(a)(dMda(a)-dnda(a)*TOF);
        
        % minimum energy ellipse data
        amin = s/2;
        
        %initial guess for a
        atilde = -amin;
        
        % Newton-Raphson solve for a
        tol = 1e-9;
        err = 1;
        while abs(err) >= tol
            a = atilde - K(atilde)/Kdot(atilde);
            err = a-atilde;
            atilde = a;
        end
        
        %calculate H2-H1
        d_anomaly = alpha(a)-beta(a);
        
        pnotrig = 4*abs(a)*(s-r1)*(s-r2)/c^2;
        pp = pnotrig*sinh((alpha(a)+beta(a))/2)^2;
        pm = pnotrig*sinh((alpha(a)-beta(a))/2)^2;
        if typeid
            p = min([pp pm]);
        else
            p = max([pp pm]);
        end
        
    case 0
        typeconic = 'P';
        % to elaborate later
        
    case 1
        % alpha equations
        alpha0 = @(a)(2*asin(sqrt(s/2/a))); %rad
        dalpha0da = @(a)(-sqrt(s/(2*a-s)/a^2));
        alpha = @(a,typeb)(2*pi*typeb+(1-2*typeb)*alpha0(a));
        dalphada = @(a,typeb)((1-2*typeb)*dalpha0da(a));
        
        % beta equations
        beta0 = @(a)(2*asin(sqrt((s-c)/2/a))); %rad
        dbeta0da = @(a)(-sqrt((s-c)/(2*a-s+c)/a^2));
        beta = @(a)(type1p*beta0(a));
        dbetada = @(a)(type1p*dbeta0da(a));
        
        % mean motion equations
        n = @(a)(sqrt(mu/a^3));
        dnda = @(a)(-sqrt(mu)*3/2*a^(-5/2));
        
        % Lambert equations
        Q = @(E)(E-sin(E));
        dQdE = @(E)(1-cos(E));
        dM = @(a,typeb)(Q(alpha(a,typeb))-Q(beta(a)));
        dMda = @(a,typeb)(dQdE(alpha(a,typeb))*dalphada(a,typeb)-...
            dQdE(beta(a))*dbetada(a));
        K = @(a,typeb)(dM(a,typeb)-n(a)*TOF);
        Kdot = @(a,typeb)(dMda(a,typeb)-dnda(a)*TOF);
        
        % minimum energy ellipse data
        amin = s/2;
        nmin = n(amin);
        dMmin = dM(amin,0); %(alpha0(amin)-sin(alpha0(amin))) - (beta0(amin)-sin(beta0(amin)));
        TOFmin = dMmin/nmin;
        
        % determine ellipse type
        switch sign(TOF - TOFmin)
            case -1
                typeconic = 'A';
                typeb = 0;
            case 0
                typeconic = 'M';
                typeb = 0;
            case 1
                typeconic = 'B';
                typeb = 1;
        end
        
        %initial guess for a
%         ntilde = dMmin/(TOFmin+abs(TOF-TOFmin));
%         atilde = ((mu/ntilde^2)^(1/3)+2*amin)/3;
        atilde = amin+1;
        
%         arange = amin:Re/10:atilde+Re;
%         for i = 1:length(arange)
%             kplot(i) = K(arange(i),typeb);
%         end
%         plot(arange,kplot,atilde,K(atilde,typeb),'*')        
        
        % Newton-Raphson solve for a
        tol = 1e-9;
        err = 1;
        while abs(err) >= tol
            a = atilde - K(atilde,typeb)/Kdot(atilde,typeb);
            err = a-atilde;
            atilde = a;
        end
        
        %calculate E2-E1
        d_anomaly = alpha(a,typeb)-beta(a); %rad
        
        pnotrig = 4*a*(s-r1)*(s-r2)/c^2;
        pp = pnotrig*sin((alpha(a,typeb)+beta(a))/2)^2;
        pm = pnotrig*sin((alpha(a,typeb)-beta(a))/2)^2;
        if xor(typeid,typeb)
            p = min([pp pm]);
        else
            p = max([pp pm]);
        end
end

type = [num2str(typenum) typeconic];

return