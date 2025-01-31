%  ------------------------------------------------------------------------------
%
%                           procedure doubler
%
%  this rountine accomplishes the iteration evaluation of one step for the double-r
%  angles only routine.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    cc1          -
%    cc2          -
%    magrsite1    -
%    magrsite2    -
%    magr1in      - initial estimate
%    magr2in      - initial estimate
%    los1         - line of sight vector for 1st
%    los2         - line of sight vector for 2nd
%    los3         - line of sight vector for 3rd
%    rsite1       - eci site position vector #1                  km
%    rsite2       - eci site position vector #2                  km
%    rsite3       - eci site position vector #3                  km
%    t1           -                                              sec
%    t3           -                                              sec
%    dm           - direction of motion                          'S', 'L'
%       likely not needed because for orbits that repeat, you always want to go the short way, rel to the initial vectors
%    n12          - # of days between tracks 1,2                 days
%    n13          - # of days between tracks 1,3                 days
%    n23          - # of days between tracks 2,3                 days
%
%  outputs        :
%    r2           - position vector at t2                        km
%    r3           - position vector at t3                        km
%    f1           -
%    f2           -
%    q1           - quality estimate
%    magr1        - magnitude of r1 vector                       km
%    magr2        - magnitude of r2 vector                       km
%    a            - semi major axis                              km / s
%    deltae32     - eccentric anomaly difference 3-2
%    errstr       - output results for debugging
%
%  coupling       :
%    dot, cross, mag
%
%  references     :
%    vallado       2022, 450 Alg 53
%
%  [r2, r3, f1, f2, q1, magr1, magr2, a, deltae32] = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,...
%        los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, n12, n13, n23);
% ----------------------------------------------------------------------------


function [r2, r3, f1, f2, q1, magr1, magr2, a, deltae32] = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,...
        los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, n12, n13, n23)

    constastro;

    %show = 'y';
    show = 'n';

    twopi = 2.0 * pi;

    %   re = 6378.137;
    %   mu = 3.986004418e5;

    %   re = 149597870.0;
    %   mu = 1.32712428e11;

    tempsq = cc1 * cc1 - 4.0 * (magrsite1^2 - magr1in^2);
    if tempsq < 0.0
        tempsq = 300.0;  % default range, use because hyperbolic likely at shorter times, lower alt
    end
    rho1 = (-cc1 + sqrt(tempsq)) * 0.5;
    tempsq = cc2 * cc2 - 4.0 * (magrsite2^2 - magr2in^2);
    if tempsq < 0.0
        tempsq = 300.0;  % a default range, not 0.0
    end
    rho2 = (-cc2 + sqrt(tempsq)) * 0.5;

    %rsite1
    %rsite2
    r1 = rho1 * los1 + rsite1;
    r2 = rho2 * los2 + rsite2;

    %rho1
    %r1

    if show == 'y'
        fprintf(1,'start of loop  %11.7f  %11.7f  \n',magr1in,magr2in );
    end
    magr1 = mag(r1);
    magr2 = mag(r2);
    if show == 'y'
        fprintf(1,'r1  %11.7f   %11.7f  %11.7f \n',r1);
        fprintf(1,'r2  %11.7f   %11.7f  %11.7f \n',r2);
    end
    %    fprintf(1,'r1  %11.7f   %11.7f  %11.7f \n',r1/re);
    %    fprintf(1,'r2  %11.7f   %11.7f  %11.7f   %11.7f  %11.7f  %11.7f \n',r2/re, r2);


    %if direct == 'y'
    w = cross(r1,r2)/(magr1*magr2);
    %else
    %    w = -cross(r1,r2)/(magr1*magr2);
    %end

    % change to negative sign
    rho3 =  -dot(rsite3,w)/dot(los3,w);
    %rho1
    %rho2
    %rho3
    %los1
    %los2
    %los3
    %pause;
    r3 = rho3.*los3 + rsite3;
    if show == 'y'
        fprintf(1,'r3  %11.7f   %11.7f  %11.7f \n',r3);
    end
    %    fprintf(1,'r3  %11.7f   %11.7f  %11.7f \n',r3/re);
    magr3 = mag(r3);
    if show == 'y'
        fprintf(1,'after 1st mag  %11.7f  %11.7f  %11.7f \n',magr1,magr2,magr3 );
    end
    %dot(r2,r1)
    %magr2*magr1
    %cross(r2,r1)
    cosdv21 = dot(r2,r1)/(magr2*magr1);
    sindv21 = mag(cross(r2,r1))/(magr2*magr1);
    dv21 = atan2(sindv21,cosdv21) + twopi * n12;

    cosdv31 = dot(r3,r1)/(magr3*magr1);
    %    sindv31 = mag(cross(r3,r1))/(magr3*magr1);
    sindv31 = sqrt(1.0 - cosdv31^2);
    dv31 = atan2(sindv31,cosdv31) + twopi * n13;

    cosdv32 = dot(r3,r2)/(magr3*magr2);
    sindv32 = mag(cross(r3,r2))/(magr3*magr2);

    dv32 = atan2(sindv32,cosdv32) + twopi * n23;
    if dv31 > pi
        c1 = (magr2*sindv32)/(magr1*sindv31);
        c3 = (magr2*sindv21)/(magr3*sindv31);
        %p = (c1*magr1+c3*magr3-magr2)/(c1+c3-1);
        p = (sindv21 - sindv31 + sindv32) / (-sindv31 / magr2 + sindv21 / magr3 + sindv32 / magr1);
    else
        c1 = (magr1*sindv31)/(magr2*sindv32);
        c3 = (magr1*sindv21)/(magr3*sindv32);
        p = (c3*magr3-c1*magr2+magr1)/(-c1+c3+1);
    end

    ecosv1 = p/magr1-1;
    ecosv2 = p/magr2-1;
    ecosv3 = p/magr3-1;

    if dv21 ~= pi
        esinv2 = (-cosdv21*ecosv2+ecosv1)/sindv21;
    else
        esinv2 = (cosdv32*ecosv2-ecosv3)/sindv32;  % 32!! not 31
    end

    e = sqrt(ecosv2^2+esinv2^2);
    a = p/(1-e^2);

    if e*e < 1.0
        n = sqrt(mu/a^3);

        s = magr2/p*sqrt(1-e^2)*esinv2;
        c = magr2/p*(e^2+ecosv2);

        sinde32 = magr3/sqrt(a*p)*sindv32-magr3/p*(1-cosdv32)*s;
        cosde32 = 1-magr2*magr3/(a*p)*(1-cosdv32);
        deltae32 = atan2(sinde32,cosde32) + twopi * n23;

        sinde21 = magr1/sqrt(a*p)*sindv21+magr1/p*(1-cosdv21)*s;
        cosde21 = 1-magr2*magr1/(a*p)*(1-cosdv21);
        deltae21 = atan2(sinde21,cosde21)+ twopi * n12;

        deltam32 = deltae32+2*s*(sin(deltae32/2))^2-c*sin(deltae32);
        deltam12 = -deltae21+2*s*(sin(deltae21/2))^2+c*sin(deltae21);
    else
        fprintf(1,'hyperbolic, e1 is greater than 0.99 %11.7f \n',e);
        if (a > 0.0)
            a = -a;  % get right sign for hyperbolic orbit
            p = -p;
        end
        n = sqrt(mu/-a^3);

        s = magr2/p*sqrt(e^2-1)*esinv2;
        c = magr2/p*(e^2+ecosv2);

        sindh32 = magr3/sqrt(-a*p)*sindv32-magr3/p*(1-cosdv32)*s;
        sindh21 = magr1/sqrt(-a*p)*sindv21+magr1/p*(1-cosdv21)*s;

        deltah32 = log( sindh32 + sqrt(sindh32^2 +1) );
        deltah21 = log( sindh21 + sqrt(sindh21^2 +1) );

        deltam32 = -deltah32+2*s*(sinh(deltah32/2))^2+c*sinh(deltah32);
        deltam12 = deltah21+2*s*(sinh(deltah21/2))^2-c*sinh(deltah21);
        % what if ends on hyperbolic solution.
        % how to pass back deltae32?
        deltae32=deltah32; %  fix
    end

    if show == 'y'
        fprintf(1,'dm32 %11.7f  dm12 %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n',deltam32,deltam12,c1,c3,p,a,e,s,c );
        fprintf(1,'%11.7f %11.7f %11.7f \n',dv21,dv31,dv32 );
    end
    f1 = t1-deltam12/n;
    f2 = t3-deltam32/n;

    q1 = sqrt(f1^2+f2^2);

end

