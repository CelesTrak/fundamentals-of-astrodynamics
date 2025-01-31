% ------------------------------------------------------------------------------
%
%                           function getGaussRoot
%
%  this function solves the problem of orbit determination using three
%    optical sightings.  the solution function uses the gaussian technique.
%    there are lots of debug statements in here to test various options.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  23 dec 2003
%   8 oct 2007
%
%  inputs          description                    range / units
%    re           - radius earth, sun, etc        km
%    mu           - grav param earth, sun etc     km3/s2
%    rtasc1       - right ascension #1                          rad
%    rtasc2       - right ascension #2                          rad
%    rtasc3       - right ascension #3                          rad
%    decl1        - declination #1                              rad
%    decl2        - declination #2                              rad
%    decl3        - declination #3                              rad
%    jd1, jdf1    - julian date of 1st sighting                 days from 4713 bc
%    jd2, jdf2    - julian date of 2nd sighting                 days from 4713 bc
%    jd3, jdf3    - julian date of 3rd sighting                 days from 4713 bc
%    rsite1       - eci site position vector                    km
%    rsite2       - eci site position vector                    km
%    rsite3       - eci site position vector                    km
%
%  outputs        :
%    r            - ijk position vector at t2     km
%    v            - ijk velocity vector at t2     km / s
%
%  locals         :
%    l1           - line of sight vector for 1st
%    l2           - line of sight vector for 2nd
%    l3           - line of sight vector for 3rd
%    tau          - taylor expansion series about
%                   tau ( t - to )
%    tausqr       - tau squared
%    t21t23       - (t2-t1) * (t2-t3)
%    t31t32       - (t3-t1) * (t3-t2)
%    i            - index
%    d            -
%    rho          - range from site to sat at t2  km
%    rhodot       -
%    dmat         -
%    rs1          - site vectors
%    rs2          -
%    rs3          -
%    earthrate    - velocity of earth rotation
%    p            -
%    q            -
%    oldr         -
%    oldv         -
%    f1           - f coefficient
%    g1           -
%    f3           -
%    g3           -
%    l2dotrs      -
%
%  coupling       :
%    mag          - magnitude of a vector
%    detrminant   - evaluate the determinant of a matrix
%    factor       - find roots of a polynomial
%    matmult      - multiply two matrices together
%    gibbs        - gibbs method of orbit determination
%    hgibbs       - herrick gibbs method of orbit determination
%    angl         - angl between two vectors
%
%  references     :
%    vallado       2007, 429-439
%
% [bigr2] = getGaussRoot(decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, ...
%                   jd1, jdf1, jd2, jdf2, jd3, jdf3, rseci1, rseci2, rseci3)
% ------------------------------------------------------------------------------

function [bigr2] = getGaussRoot(decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, ...
        jd1, jdf1, jd2, jdf2, jd3, jdf3, rseci1, rseci2, rseci3)
    show = 'n';
    constastro;
    rad = 180.0 / pi;

    % -------------------------  implementation   -------------------------
    tu = sqrt(re * re * re / mu)

    % ---------- set middle to 0, find decls to others -----------
    tau12 = (jd1 - jd2) * 86400.0 + (jdf1 - jdf2) * 86400.0; % days to sec
    tau13 = (jd1 - jd3) * 86400.0 + (jdf1 - jdf3) * 86400.0;
    tau32 = (jd3 - jd2) * 86400.0 + (jdf3 - jdf2) * 86400.0;

    tau12c = tau12 / tu;  % to tu
    tau13c = tau13 / tu;
    tau32c = tau32 / tu;

    fprintf(1,'jd123 %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f  \n',jd1,jdf1, jd2,jdf2, jd3, jdf3);
    fprintf(1,'tau12 %14.6f tau13  %14.6f tau32  %14.6f \n',tau12,tau13, tau32);

    % ----------------  find line of sight unit vectors  ---------------
    los1(1)= cos(decl1)*cos(rtasc1);
    los1(2)= cos(decl1)*sin(rtasc1);
    los1(3)= sin(decl1);

    los2(1)= cos(decl2)*cos(rtasc2);
    los2(2)= cos(decl2)*sin(rtasc2);
    los2(3)= sin(decl2);

    los3(1)= cos(decl3)*cos(rtasc3);
    los3(2)= cos(decl3)*sin(rtasc3);
    los3(3)= sin(decl3);

    % ------------- find l matrix and determinant -----------------
    if show == 'y'
        los1
    end
    vs = [0 0 0]';
    aecef = [0 0 0]';
    %[l1eci,vs3,aeci] = ecef2eci(l1',vs,aecef,(jd1-2451545.0)/36525.0,jd1,0.0,0.0,0.0,0,ddpsi,ddeps);
    %[l2eci,vs3,aeci] = ecef2eci(l2',vs,aecef,(jd2-2451545.0)/36525.0,jd2,0.0,0.0,0.0,0,ddpsi,ddeps);
    %[l3eci,vs3,aeci] = ecef2eci(l3',vs,aecef,(jd3-2451545.0)/36525.0,jd3,0.0,0.0,0.0,0,ddpsi,ddeps);

    l1eci = los1;
    l2eci = los2;
    l3eci = los3;
    % leave these as they come since the topoc radec are already eci
    if show == 'y'
        l1eci
    end
    % --------- called lmati since it is only used for determ -----
    for i= 1 : 3
        lmat(i,1) = l1eci(i);
        lmat(i,2) = l2eci(i);
        lmat(i,3) = l3eci(i);
        rsmatc(i,1)  = rseci1(i)/re;  % eci in er
        rsmatc(i,2)  = rseci2(i)/re;
        rsmatc(i,3)  = rseci3(i)/re;
        rs2c(i) = rseci2(i) / re;  % er for later
    end
    if show == 'y'
        lmat

        fprintf(1,'rsmat eci %11.7f  %11.7f  %11.7f km \n',rsmatc');
    end

    fprintf(1,'this should be the inverse of what the code finds later\n');
    li = inv(lmat)


    %     % alt way of Curtis not seem to work ------------------
    %     p1 = cross(los2, los3);
    %     p2 = cross(los1, los3);
    %     p3 = cross(los1, los2);
    %     % both are the same
    %     dx = dot(los1,p1);
    %     %dx = dot(los3,p3);
    %     lcmat(1,1) = dot(rsite1,p1);
    %     lcmat(2,1) = dot(rsite2,p1);
    %     lcmat(3,1) = dot(rsite3,p1);
    %     lcmat(1,2) = dot(rsite1,p2);
    %     lcmat(2,2) = dot(rsite2,p2);
    %     lcmat(3,2) = dot(rsite3,p2);
    %     lcmat(1,3) = dot(rsite1,p3);
    %     lcmat(2,3) = dot(rsite2,p3);
    %     lcmat(3,3) = dot(rsite3,p3);
    %     tau31= (jd3-jd1)*86400.0;
    %     if show == 'y'
    %         lcmat
    %     end
    %     aa = 1/dx *(-lcmat(1,2)*tau32c/tau13c + lcmat(2,2) + lcmat(3,2)*tau12c/tau13c);
    %     bb = 1/(6.0*dx) *(lcmat(1,2)*(tau32c*tau32c - tau13c*tau13c)*tau32c/tau13c  ...
    %         + lcmat(3,2)*(tau13c*tau13c - tau12c*tau12c) * tau12c/tau13c);
    %     % alt way of Curtis not seem to work ------------------


    d= det(lmat);
    if show == 'y'
        d
    end
    % ------------------ now assign the inverse -------------------
    lmati(1,1) = ( l2eci(2)*l3eci(3)-l2eci(3)*l3eci(2)) / d;
    lmati(2,1) = (-l1eci(2)*l3eci(3)+l1eci(3)*l3eci(2)) / d;
    lmati(3,1) = ( l1eci(2)*l2eci(3)-l1eci(3)*l2eci(2)) / d;
    lmati(1,2) = (-l2eci(1)*l3eci(3)+l2eci(3)*l3eci(1)) / d;
    lmati(2,2) = ( l1eci(1)*l3eci(3)-l1eci(3)*l3eci(1)) / d;
    lmati(3,2) = (-l1eci(1)*l2eci(3)+l1eci(3)*l2eci(1)) / d;
    lmati(1,3) = ( l2eci(1)*l3eci(2)-l2eci(2)*l3eci(1)) / d;
    lmati(2,3) = (-l1eci(1)*l3eci(2)+l1eci(2)*l3eci(1)) / d;
    lmati(3,3) = ( l1eci(1)*l2eci(2)-l1eci(2)*l2eci(1)) / d;
    if show == 'y'
        lmati
    end

    lir = lmati*rsmatc;

    % ------------ find f and g series at 1st and 3rd obs ---------
    %   speed by assuming circ sat vel for udot here ??
    %   some similartities in 1/6t3t1 ...
    % --- keep separated this time ----
    a1c =  tau32c / (tau32c - tau12c);
    a1uc=  (tau32c*((tau32c-tau12c)^2 - tau32c*tau32c )) / (6.0*(tau32c - tau12c));
    a3c = -tau12c / (tau32c - tau12c);
    a3uc= -(tau12c*((tau32c-tau12c)^2 - tau12c*tau12c )) / (6.0*(tau32c - tau12c));

    if show == 'y'
        fprintf(1,'a1/a3 %11.7f  %11.7f  %11.7f  %11.7f \n',a1c,a1uc,a3c,a3uc );
    end
    % --- form initial guess of r2 ----
    d1c=  lir(2,1)*a1c - lir(2,2) + lir(2,3)*a3c;
    d2c=  lir(2,1)*a1uc + lir(2,3)*a3uc;
    if show == 'y'
        lir
        d1c
        d2c
    end

    % ------- solve eighth order poly not same as laplace ---------
    magrs2c = mag(rseci2) / re;
    l2dotrs= dot( los2,rseci2 ) / re;
    if show == 'y'
        fprintf(1,'magrs2 %11.7f  %11.7f  \n',magrs2c,l2dotrs );
    end
    muc = 1.0;

    poly( 1)=  1.0;  % r2^8th variable%%%%%%%%%%%%%%
    poly( 2)=  0.0;
    poly( 3)=  -(d1c*d1c + 2.0*d1c*l2dotrs + magrs2c^2);
    poly( 4)=  0.0;
    poly( 5)=  0.0;
    poly( 6)=  -2.0*muc*(l2dotrs*d2c + d1c*d2c);
    poly( 7)=  0.0;
    poly( 8)=  0.0;
    poly( 9)=  -muc*muc*d2c*d2c;
    if show == 'y'
        fprintf(1,'%11.7f \n',poly);
    end

    rootarr = roots( poly );
    if show == 'y'
        rootarr;
        %fprintf(1,'rootarr %11.7f \n',rootarr);
    end

    % ------------------ select the correct root ------------------
    bigr2c= -99999990.0;
    % change from 1
    for j= 1 : 8
        if ( rootarr(j) > bigr2c ) & ( isreal(rootarr(j)) )
            bigr2c= rootarr(j);
        end  % if (
    end

    % simply iterate to find the correct root
    bigr2 = 100.0;
    % can do Newton iteration Curtis p 289 simple derivative of the poly
    % makes sense since the values are so huge in the polynomial
    % do as Halley iteration since derivatives are possible. tests at LEO, GPS, GEO,
    % all seem to converge to the proper answer
    kk = 0;
    bigr2c = 20000.0 / re; % er guess ~GPS altitude
    while (abs(bigr2 - bigr2c) > 0.5 && kk < 15)
        bigr2 = bigr2c;
        deriv = bigr2^8 + poly(3) * bigr2^6 + poly(6) * bigr2^3 + poly(9);
        derivn1 = 8.0 * bigr2^7 + 6.0 * poly(3) * bigr2^5 + 3.0 * poly(6) * bigr2^2;
        derivn2 = 56.0 * bigr2^6 + 30.0 * poly(3) * bigr2^4 + 6.0 * poly(6) * bigr2;
        % Newton iteration
        %bigr2n = bigr2 - deriv / derivn1;
        % Halley iteration
        bigr2c = bigr2 - (2.0 * deriv * derivn1) / (2.0 * derivn1 * derivn1 - deriv * derivn2);
        kk = kk + 1;
    end

    bigr2 = bigr2c * re;
    if show == 'y'
        bigr2;
    end

end



