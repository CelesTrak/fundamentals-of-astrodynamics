%  ------------------------------------------------------------------------------
%
%                           procedure anglesgauss
%
%  this procedure solves the problem of orbit determination using three
%    optical sightings. the solution procedure uses the gaussian technique.
%    the 8th order root is generally the big point of discussion. A Halley iteration
%    permits a quick solution to find the correct root, with a starting guess of 20000 km.
%    the general formulation yields polynomial coefficients that are very large, and can easily
%    become overflow operations. Thus, canonical units are used only until the root is found,
%    then regular units are resumed.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    tdecl1       - declination #1                               rad
%    tdecl2       - declination #2                               rad
%    tdecl3       - declination #3                               rad
%    trtasc1      - right ascension #1                           rad
%    trtasc2      - right ascension #2                           rad
%    trtasc3      - right ascension #3                           rad
%    jd1, jdf1    - julian date of 1st sighting                  days from 4713 bc
%    jd2, jdf2    - julian date of 2nd sighting                  days from 4713 bc
%    jd3, jdf3    - julian date of 3rd sighting                  days from 4713 bc
%    rs1          - eci site position vector #1                  km
%    rs2          - eci site position vector #2                  km
%    rs3          - eci site position vector #3                  km
%
%  outputs        :
%    r2           -  position vector at t2                       km
%    v2           -  velocity vector at t2                       km / s
%    errstr       - output results for debugging
%
%  locals         :
%    los1         - line of sight vector for 1st
%    los2         - line of sight vector for 2nd
%    los3         - line of sight vector for 3rd
%    tau          - taylor expansion series about
%                   tau ( t - to )
%    tausqr       - tau squared
%    t21t23       - (t2-t1) * (t2-t3)
%    t31t32       - (t3-t1) * (t3-t2)
%    i            - index
%    d            -
%    rho          - range from site to sat at t2                 km
%    rhodot       -
%    dmat         -
%    earthrate    - velocity of earth rotation
%
%  coupling       :
%    mag          - magnitude of a vector
%    determinant  - evaluate the determinant of a matrix
%    factor       - find roots of a polynomial
%    matmult      - multiply two matrices together
%    gibbs        - gibbs method of orbit determination
%    hgibbs       - herrick gibbs method of orbit determination
%    angle        - angle between two vectors
%
%  references     :
%    vallado       2022, 448, alg 52, ex 7-2
%
% [r2, v2] = anglesg(decl1, decl2, decl3, rtasc1, rtasc2, ...
%        rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, diffsites, rseci1, rseci2, rseci3, fid);
% ------------------------------------------------------------------------------

function [r2, v2] = anglesg(decl1, decl2, decl3, rtasc1, rtasc2, ...
        rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, diffsites, rseci1, rseci2, rseci3, fid)
    %show = 'y';
    show = 'n';
    constastro;
    rad = 180.0 / pi;

    % -------------------------  implementation   -------------------------
    ddpsi = 0.0;  % delta corrections not needed for this level of precision
    ddeps = 0.0;
    magr1in = 2.0*re; % initial guesses
    magr2in = 2.01*re;
    direct = 'y';  % direction of motion short way

    % ---------- set middle to 0, find decls to others -----------
    tau12 = (jd1 - jd2) * 86400.0 + (jdf1 - jdf2) * 86400.0; % days to sec
    tau13 = (jd1 - jd3) * 86400.0 + (jdf1 - jdf3) * 86400.0;
    tau32 = (jd3 - jd2) * 86400.0 + (jdf3 - jdf2) * 86400.0;
    fprintf(1,'jd123 %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f  \n',jd1,jdf1, jd2,jdf2, jd3, jdf3);
    fprintf(1,'tau12 %14.6f tau13  %14.6f tau32  %14.6f \n',tau12,tau13, tau32);

    % switch to cannonical for better accuracy in roots
    tau12c = tau12 / tusec;
    tau13c = tau13 / tusec;
    tau32c = tau32 / tusec;

    % switch to canonical
    for i=1:3
        rs1c(i) = rseci1(i) / re;
        rs2c(i) = rseci2(i) / re;
        rs3c(i) = rseci3(i) / re;
    end

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

    % topo to body fixed (ecef)
    %     latgd = 40.0/rad;
    %     lon   = -110.0/rad;
    %     l1
    %     outv = rot2(l1, -pi + latgd);
    %     l1 = rot3(outv, -lon);
    %     outv = rot2(l2, -pi + latgd);
    %     l2 = rot3(outv, -lon);
    %     outv = rot2(l3, -pi + latgd);
    %     l3 = rot3(outv, -lon);
    %     l1
    %
    %     % take the middle trans from eecef to eci
    %     tm =  [-0.830668624503591  -0.556765707115059   0.001258429966118;...
    %    0.556766123167794  -0.830669298539751  -0.000023583565998;...
    %    0.001058469658016   0.000681061045186   0.999999207898604];
    %
    %
    % l1 = tm*l1';
    % l2=tm*l2';
    % l3=tm*l3';

    % ------------- find l matrix and determinant -----------------
    if show == 'y'
        los1
    end
    vs = [0 0 0]';
    aecef = [0 0 0]';
    %[los1,vs3,aeci] = ecef2eci(l1',vs,aecef,(jd1-2451545.0)/36525.0,jd1,0.0,0.0,0.0,0,ddpsi,ddeps);
    %[los2,vs3,aeci] = ecef2eci(l2',vs,aecef,(jd2-2451545.0)/36525.0,jd2,0.0,0.0,0.0,0,ddpsi,ddeps);
    %[los3,vs3,aeci] = ecef2eci(l3',vs,aecef,(jd3-2451545.0)/36525.0,jd3,0.0,0.0,0.0,0,ddpsi,ddeps);

    % leave these as they come since the topoc radec are already eci
    if show == 'y'
        los1
    end
    % --------- called lmati since it is only used for determ -----
    for i= 1 : 3
        lmat(i,1) = los1(i);
        lmat(i,2) = los2(i);
        lmat(i,3) = los3(i);

        rsmat(i,1)  = rseci1(i);  % eci
        rsmat(i,2)  = rseci2(i);
        rsmat(i,3)  = rseci3(i);

        rsmatc(i,1)  = rseci1(i)/re;  % eci
        rsmatc(i,2)  = rseci2(i)/re;
        rsmatc(i,3)  = rseci3(i)/re;
        rs2c(i) = rseci2(i) / re;  % er for later
    end
    if show == 'y'
        lmat

        fprintf(1,'rsmat eci %11.7f  %11.7f  %11.7f km \n',rsmat');

        % the order is right, but to print out, need '
        fprintf(1,'rsmat eci %11.7f  %11.7f  %11.7f \n',rsmat'/re);

        lmat
        fprintf(1,'this should be the inverse of what the code finds later\n');
        li = inv(lmat)
    end
    % alt way of Curtis not seem to work ------------------
    % p1 = cross(los2, los3);
    % p2 = cross(los1, los3);
    % p3 = cross(los1, los2);
    % % both are the same
    % dx = dot(los1,p1);
    % %dx = dot(los3,p3);
    % lcmat(1,1) = dot(rseci1,p1);
    % lcmat(2,1) = dot(rseci2,p1);
    % lcmat(3,1) = dot(rseci3,p1);
    % lcmat(1,2) = dot(rseci1,p2);
    % lcmat(2,2) = dot(rseci2,p2);
    % lcmat(3,2) = dot(rseci3,p2);
    % lcmat(1,3) = dot(rseci1,p3);
    % lcmat(2,3) = dot(rseci2,p3);
    % lcmat(3,3) = dot(rseci3,p3);
    % tau31= (jd3-jd1)*86400.0;
    % if show == 'y'
    %     lcmat
    % end
    % aa = 1/dx *(-lcmat(1,2)*tau32/tau31 + lcmat(2,2) + lcmat(3,2)*tau12/tau31)
    % bb = 1/(6.0*dx) *(lcmat(1,2)*(tau32*tau32 - tau31*tau31)*tau32/tau31  ...
    %     + lcmat(3,2)*(tau31*tau31 - tau12*tau12) * tau12/tau31)
    % alt way of Curtis not seem to work ------------------


    d= det(lmat);
    if show == 'y'
        d
    end
    % ------------------ now assign the inverse -------------------
    lmati(1,1) = ( los2(2)*los3(3)-los2(3)*los3(2)) / d;
    lmati(2,1) = (-los1(2)*los3(3)+los1(3)*los3(2)) / d;
    lmati(3,1) = ( los1(2)*los2(3)-los1(3)*los2(2)) / d;
    lmati(1,2) = (-los2(1)*los3(3)+los2(3)*los3(1)) / d;
    lmati(2,2) = ( los1(1)*los3(3)-los1(3)*los3(1)) / d;
    lmati(3,2) = (-los1(1)*los2(3)+los1(3)*los2(1)) / d;
    lmati(1,3) = ( los2(1)*los3(2)-los2(2)*los3(1)) / d;
    lmati(2,3) = (-los1(1)*los3(2)+los1(2)*los3(1)) / d;
    lmati(3,3) = ( los1(1)*los2(2)-los1(2)*los2(1)) / d;
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
        dl1
        dl2
    end

    % ------- solve eighth order poly not same as laplace ---------
    % switch to canonical to prevent overflows in the poly
    magrs2 = mag(rs2c);
    l2dotrs= dot( los2,rs2c );
    if show == 'y'
        fprintf(1,'magrs2 %11.7f  %11.7f  \n',magrs2,l2dotrs );
    end

    poly( 1)=  1.0;  % r2^8th variable%%%%%%%%%%%%%%
    poly( 2)=  0.0;
    poly( 3)=  -(d1c*d1c + 2.0*d1c*l2dotrs + magrs2^2);
    poly( 4)=  0.0;
    poly( 5)=  0.0;
    poly( 6)=  -2.0*(l2dotrs*d2c + d1c*d2c);
    poly( 7)=  0.0;
    poly( 8)=  0.0;
    poly( 9)=  -d2c*d2c;
    fprintf(1,'%11.7f \n',poly);

    rootarr = roots( poly );
    if show == 'y'
        rootarr;
        %fprintf(1,'rootarr %11.7f \n',rootarr);
    end

    % % ------------------ select the correct root ------------------
    % bigr2= -99999990.0;
    % % change from 1
    % for j= 1 : 8
    %     if ( rootarr(j) > bigr2 ) & ( isreal(rootarr(j)) )
    %         bigr2= rootarr(j);
    %     end  % if (
    % end
    % if show == 'y'
    %     bigr2
    % end

    bigr2 = 100.0;
    % can do Newton iteration Curtis p 289 simple derivative of the poly
    % makes sense since the values are so huge in the polynomial
    % do as Halley iteration since derivatives are possible. tests at LEO, GPS, GEO,
    % all seem to converge to the proper answer
    kk = 1;
    bigr2c = 20000.0 / re; % er guess ~GPS altitude
    % bigr2nx = bigr2c;
    while (abs(bigr2 - bigr2c) > 8.0e-5 && kk < 15)  % do in er, 0.5 km

        bigr2 = bigr2c;
        %    bigr2x = bigr2nx;
        deriv = bigr2^8 + poly(3) * bigr2^6 + poly(6) * bigr2^3 + poly(9);
        deriv1 = 8.0 * bigr2^7 + 6.0 * poly(3) * bigr2^5 + 3.0 * poly(6) *bigr2^2;
        deriv2 = 56.0 * bigr2^6 + 30.0 * poly(3) * bigr2^4 + 6.0 * poly(6) * bigr2;
        % just use Halley
        %bigr2n = bigr2 - deriv * n2 / (deriv1 * (n2 - deriv * deriv2 * 0.5));
        % Halley iteration
        bigr2c = bigr2 - (2.0 * deriv * deriv1) / (2.0 * deriv1 * deriv1 - deriv * deriv2);
    end

    if (bigr2c < 0.0 || bigr2c * re > 50000.0)
        bigr2c = 35000.0 / re;  % simply set this to about GEO, allowing for less than that too.
    end
    % now back to normal units
    bigr2 = bigr2c * re;  % km
    %    bigr2nx = bigr2nx * gravConst.re;

    a1 = tau32 / (tau32 - tau12);
    a1u = (tau32 * ((tau32 - tau12) * (tau32 - tau12) - tau32 * tau32)) /...
        (6.0 * (tau32 - tau12));
    a3 = -tau12 / (tau32 - tau12);
    a3u = -(tau12 * ((tau32 - tau12) * (tau32 - tau12) - tau12 * tau12)) /...
        (6.0 * (tau32 - tau12));
    lir = lmati * rsmat;


    % ------------ solve matrix with u2 better known --------------
    u= mu / ( bigr2*bigr2*bigr2 );

    c1= a1 + a1u*u;
    c2 = -1.0;
    c3= a3 + a3u*u;

    if show == 'y'
        fprintf(1,'u %17.14f c1 %11.7f c3 %11.7f %11.7f \n',u,c1,c2,c3);
    end

    cmat(1,1)= -c1;
    cmat(2,1)= -c2;
    cmat(3,1)= -c3;
    rhomat = lir*cmat;

    rhoold1=  rhomat(1,1)/c1;
    rhoold2=  rhomat(2,1)/c2;
    rhoold3=  rhomat(3,1)/c3;
    if show == 'y'
        fprintf(1,'rhoold %11.7f %11.7f %11.7f \n',rhoold1,rhoold2,rhoold3);
        %   fprintf(1,'rhoold %11.7f %11.7f %11.7f \n',rhoold1/re,rhoold2/re,rhoold3/re);
    end

    % ---- now form the three position vectors -----
    for i= 1 : 3
        r1(i)=  rhomat(1,1)*los1(i)/c1 + rseci1(i);
        r2(i)=  rhomat(2,1)*los2(i)/c2 + rseci2(i);
        r3(i)=  rhomat(3,1)*los3(i)/c3 + rseci3(i);
    end
    if show == 'y'
        fprintf(1,'r1 %11.7f %11.7f %11.7f \n',r1);
        fprintf(1,'r2 %11.7f %11.7f %11.7f \n',r2);
        fprintf(1,'r3 %11.7f %11.7f %11.7f \n',r3);
    end

    % now find the middle velocity vector with gibbs or hgibbs from end of formal Gauss
    [v2, theta, theta1, copa, error] = gibbs(r1, r2, r3);

    %pause;
    % -------- loop through the refining process ------------  while () for
    if show == 'y'
        fprintf(1,'now refine the answer \n');
    end
    rho2 = 999999.9;
    ll   = 0;
    % disabled now...
    while ((abs(rhoold2-rho2) > 1.0e-12) && (ll <= -1 ))  % ll <= 4   15
        ll = ll + 1;
        fprintf(1, ' iteration #%3i \n',ll );
        rho2 = rhoold2;  % reset now that inside while loop

        % ---------- now form the three position vectors ----------
        for i= 1 : 3
            r1(i)=  rhomat(1,1)*los1(i)/c1 + rseci1(i);
            r2(i)= -rhomat(2,1)*los2(i)    + rseci2(i);
            r3(i)=  rhomat(3,1)*los3(i)/c3 + rseci3(i);
        end
        magr1 = mag( r1 );
        magr2 = mag( r2 );
        magr3 = mag( r3 );

        [v2,theta,theta1,copa,error] = gibbsh(r1,r2,r3, re, mu);

        rad = 180.0/pi;
        fprintf(1,'r1 %16.14f %16.14f %16.14f %11.7f %11.7f %16.14f \n',r1,theta*rad,theta1*rad, copa*rad);
        fprintf(1,'r2 %11.7f %11.7f %11.7f \n',r2);
        fprintf(1,'r3 %11.7f %11.7f %11.7f \n',r3);
        fprintf(1,'w gibbs km/s       v2 %11.7f %11.7f %11.7f \n',v2);

        % check if too close obs
        if ( (strcmp(error, '          ok') == 0) && ((abs(theta) < 1.0/rad) || (abs(theta1) < 1.0/rad)) ) % 0 is false
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2,v2, re, mu);
            if show == 'y'
                fprintf(1,'coes init ans %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f\n',...
                    p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
            end
            % --- hgibbs to get middle vector ----
            [v2,theta,theta1,copa,error] = hgibbs(r1,r2,r3,jd1+jdf1,jd2+jdf2,jd3+jdf3);
            if show == 'y'
                fprintf(1,'using hgibbs: ' );
            end
        end

        [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2,v2, re, mu);
        if show == 'y'
            fprintf(1,'coes init ans %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f\n',...
                p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
            %fprintf(1,'dr %11.7f m %11.7f m/s \n',1000*mag(r2-r2ans),1000*mag(v2-v2ans) );
        end

        if ( ll <= 8 )  % 4
            % --- now get an improved estimate of the f and g series --
            u= mu / ( magr2*magr2*magr2 );
            rdot= dot(r2,v2)/magr2;
            udot= (-3.0*mu*rdot) / (magr2^4);

            if show == 'y'
                fprintf(1,'u %17.15f rdot  %11.7f udot %11.7f \n',u,rdot,udot );
            end
            tausqr= tau12*tau12;
            f1=  1.0 - 0.5*u*tausqr -(1.0/6.0)*udot*tausqr*tau12;
            %                       - (1.0/24.0) * u*u*tausqr*tausqr
            %                       - (1.0/30.0)*u*udot*tausqr*tausqr*tau1;
            g1= tau12 - (1.0/6.0)*u*tau12*tausqr - (1.0/12.0) * udot*tausqr*tausqr;
            %                       - (1.0/120.0)*u*u*tausqr*tausqr*tau1
            %                       - (1.0/120.0)*u*udot*tausqr*tausqr*tausqr;
            tausqr= tau32*tau32;
            f3=  1.0 - 0.5*u*tausqr -(1.0/6.0)*udot*tausqr*tau32;
            %                       - (1.0/24.0) * u*u*tausqr*tausqr
            %                       - (1.0/30.0)*u*udot*tausqr*tausqr*tau3;
            g3= tau32 - (1.0/6.0)*u*tau32*tausqr - (1.0/12.0) * udot*tausqr*tausqr;
            %                       - (1.0/120.0)*u*u*tausqr*tausqr*tau3
            %                       - (1.0/120.0)*u*udot*tausqr*tausqr*tausqr;
            if show == 'y'
                fprintf(1,'f1 %11.7f g1 %11.7f f3 %11.7f g3 %11.7f \n',f1,g1,f3,g3 );
            end
        else
            % -------- use exact method to find f and g -----------
            theta  = angl( r1,r2 );
            theta1 = angl( r2,r3 );

            f1= 1.0 - ( (magr1*(1.0 - cos(theta)) / p ) );
            g1= ( magr1*magr2*sin(-theta) ) / sqrt( p );  % - angl because backwards
            f3= 1.0 - ( (magr3*(1.0 - cos(theta1)) / p ) );
            g3= ( magr3*magr2*sin(theta1) ) / sqrt( p );

        end

        c1=  g3 / (f1*g3 - f3*g1);
        c3= -g1 / (f1*g3 - f3*g1);

        if show == 'y'
            fprintf(1,' c1 %11.7f c3 %11.7f %11.7f \n',c1,c2,c3);
        end

        % ----- solve for all three ranges via matrix equation ----
        cmat(1,1)= -c1;
        cmat(2,1)= -c2;
        cmat(3,1)= -c3;
        rhomat = lir*cmat;

        if show == 'y'
            fprintf(1,'rhomat %11.7f %11.7f %11.7f \n',rhomat);
            %        fprintf(1,'rhomat %11.7f %11.7f %11.7f \n',rhomat/re);
        end

        rhoold1=  rhomat(1,1)/c1;
        rhoold2=  rhomat(2,1)/c2;
        rhoold3=  rhomat(3,1)/c3;
        if show == 'y'
            fprintf(1,'rhoold %11.7f %11.7f %11.7f \n',rhoold1,rhoold2,rhoold3);
            %   fprintf(1,'rhoold %11.7f %11.7f %11.7f \n',rhoold1/re,rhoold2/re,rhoold3/re);
        end

        for i= 1 : 3
            r1(i)=  rhomat(1,1)*los1(i)/c1 + rseci1(i);
            r2(i)=  rhomat(2,1)*los2(i)/c2 + rseci2(i);
            r3(i)=  rhomat(3,1)*los3(i)/c3 + rseci3(i);
        end
        if show == 'y'
            fprintf(1,'r1 %11.7f %11.7f %11.7f \n',r1);
            fprintf(1,'r2 %11.7f %11.7f %11.7f \n',r2);
            fprintf(1,'r3 %11.7f %11.7f %11.7f \n',r3);

            fprintf(1,'====================next loop \n');
        end
        % ----------------- check for convergence -----------------
        %pause
        if show == 'y'

            fprintf(1,'rhoold while  %16.14f %16.14f \n',rhoold2,rho2);
        end
    end   % while the ranges are still changing

    % ---------------- find all three vectors ri ------------------
    for i= 1 : 3
        r1(i)=  rhomat(1,1)*los1(i)/c1 + rseci1(i);
        r2(i)= -rhomat(2,1)*los2(i)    + rseci2(i);
        r3(i)=  rhomat(3,1)*los3(i)/c3 + rseci3(i);
    end

end