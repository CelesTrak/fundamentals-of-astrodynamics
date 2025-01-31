%
%
%
%
%
function [pnm, ppnm, accel] = learnorm(mu, rbar, recef, c, s, nmax, mmax)
    for n = 2:nmax %RAE
        norm1(n) = sqrt((2*n+1)/(2*n-1)); %RAE
        norm2(n) = sqrt((2*n+1)/(2*n-3)); %RAE
        norm11(n) = sqrt((2*n+1)/(2*n))/(2*n-1); %RAE
        for m = 1:n %RAE
            norm1m(n,m) = sqrt((n-m)*(2*n+1)/((n+m)*(2*n-1))); %RAE
            norm2m(n,m) = sqrt((n-m)*(n-m-1)*(2*n+1)/((n+m)*(n+m-1)*(2*n-3))); %RAE
        end %RAE
    end %RAE
    for n=2:nmax %RAE
        pnm(n-1,n)=0;
    end
    %recef=rnp*recef; %RAE
    e1=recef(1)^2+recef(2)^2;
    r2=e1+recef(3)^2;
    absr=sqrt(r2);
    r1=sqrt(e1);
    sphi=recef(3)/absr;
    cphi=r1/absr;
    sm(1)=0;
    cm(1)=1;
    if (r1~=0)
        sm(1)=recef(2)/r1;

        cm(1)=recef(1)/r1;
    end
    rb(1)=rbar/absr;
    rb(2)=rb(1)^2;
    sm(2)=2*cm(1)*sm(1);
    cm(2)=2*cm(1)^2-1;
    root3=sqrt(3); %RAE
    root5=sqrt(5); %RAE
    pn(1)=root3*sphi; %norm
    pn(2)=root5*(3*sphi^2-1)/2; %norm
    ppn(1)=root3; %norm
    ppn(2)=root5*3*sphi; %norm
    pnm(1,1)=root3; %norm
    pnm(2,2)=root5*root3*cphi/2; %norm
    pnm(2,1)=root5*root3*sphi; %norm %RAE
    ppnm(1,1)=-root3*sphi; %norm
    ppnm(2,2)=-root3*root5*sphi*cphi; %norm
    ppnm(2,1)=root5*root3*(1-2*sphi^2); %norm
    if (nmax>=3) %RAE
        for n=3:nmax %RAE
            nm1=n-1;
            nm2=n-2;
            rb(n)=rb(nm1)*rb(1);
            sm(n)=2*cm(1)*sm(nm1)-sm(nm2);
            cm(n)=2*cm(1)*cm(nm1)-cm(nm2);
            e1=2*n-1;
            pn(n)=(e1*sphi*norm1(n)*pn(nm1)-nm1*norm2(n)*pn(nm2))/n; %norm
            ppn(n)=norm1(n)*(sphi*ppn(nm1)+n*pn(nm1)); %norm
            pnm(n,n)=e1*cphi*norm11(n)*pnm(nm1,nm1); %norm
            ppnm(n,n)=-n*sphi*pnm(n,n);
        end
        for n=3:nmax %RAE
            nm1=n-1;
            e1=(2*n-1)*sphi;
            e2=-n*sphi;
            for m=1:nm1
                e3=norm1m(n,m)*pnm(nm1,m); %norm
                e4=n+m;
                e5=(e1*e3-(e4-1)*norm2m(n,m)*pnm(n-2,m))/(n-m); %norm
                pnm(n,m)=e5;
                ppnm(n,m)=e2*e5+e4*e3;
            end
        end
    end
    asph(1)=-1;
    asph(3)=0;
    for n=2:nmax %RAE
        ni=n+1; %RAE
        e1=c(ni,1)*rb(n); %RAE
        asph(1)=asph(1)-(n+1)*e1*pn(n);
        asph(3)=asph(3)+e1*ppn(n);
    end
    asph(3)=cphi*asph(3);
    t1=0;
    t3=0;
    asph(2)=0;
    for n=2:nmax %RAE
        ni=n+1; %RAE
        e1=0;
        e2=0;
        e3=0;
        if (n>mmax) %RAE
            nmodel=mmax; %RAE
        else %RAE
            nmodel=n; %RAE
        end %RAE
        for m=1:nmodel %RAE
            mi=m+1; %RAE
            tsnm=s(ni,mi); %RAE
            tcnm=c(ni,mi); %RAE
            tsm=sm(m);
            tcm=cm(m);
            tpnm=pnm(n,m);
            e4=tsnm*tsm+tcnm*tcm;
            e1=e1+e4*tpnm;
            e2=e2+m*(tsnm*tcm-tcnm*tsm)*tpnm;
            e3=e3+e4*ppnm(n,m);
        end
        t1=t1+(n+1)*rb(n)*e1;
        asph(2)=asph(2)+rb(n)*e2;
        t3=t3+rb(n)*e3;
    end
    e4=mu/r2;
    asph(1)=e4*(asph(1)-cphi*t1);
    asph(2)=e4*asph(2);
    asph(3)=e4*(asph(3)+t3);
    e5=asph(1)*cphi-asph(3)*sphi;
    accel(1,1)=e5*cm(1)-asph(2)*sm(1); %RAE
    accel(2,1)=e5*sm(1)+asph(2)*cm(1); %RAE
    accel(3,1)=asph(1)*sphi+asph(3)*cphi; %RAE
    %accel=rnp'*agr; %RAE
    return

