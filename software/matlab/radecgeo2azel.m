% radecgeo2azel
%
% this function finds the az-el values given geocentric rtasc decl
% for scanning geo orbits
%
%  [az, el] = radecgeo2azel(rtasc, decl, latgd, lon, alt, iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps)
%

function [az, el] = radecgeo2azel(rtasc, decl, latgd, lon, alt, iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps)
    constastro;
    rad = 180.0/pi;

    % find slant range distance to GEO sat

    % ----------------- get site vector in ecef -------------------
    [rs, vs] = site ( latgd, lon, alt );
    % find GEO location that is being looked at
    rr = 42164.0;  % km
    %     ref(1)= rr*cos(0.0)*cos(rtasc);
    %     ref(2)= rr*cos(0.0)*sin(rtasc);
    %     ref(3)= rr*sin(0.0);
    %
    %     rr = mag(rs-ref);

    %    [r,v] = radec2rv( rr,rtasc,decl,drr,drtasc,ddecl );
    reci(1)= rr*cos(decl)*cos(rtasc);
    reci(2)= rr*cos(decl)*sin(rtasc);
    reci(3)= rr*sin(decl);
    reci = reci';


    %    [rho,az,el,drho,daz,del] = rv2razel ( reci,veci, latgd,lon,alt,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );
    halfpi = pi*0.5;
    small  = 0.00000001;
    % --------------------- implementation ------------------------
    % -------------------- convert eci to ecef --------------------
    %    a = [0;0;0];
    %    [recef,vecef,aecef] = eci2ecef(reci,veci,a,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps);
    [fArgs] = fundarg(ttt, '80');

    [prec,psia,wa,ea,xa] = precess ( ttt, '80' );

    [deltapsi, trueeps, meaneps, nut] = nutation  (ttt, ddpsi, ddeps, iau80arr, fArgs);

    [st, stdot] = sidereal(jdut1, deltapsi, meaneps, fArgs(5), lod, 2, '80' );

    [pm] = polarm(xp,yp,ttt,'80');

    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    rpef  = st'*nut'*prec'*reci;
    recef = pm'*rpef;

    % ------- find ecef range vector from site to satellite -------
    rhoecef  = recef - rs;
    %    drhoecef = vecef;
    %    rho      = mag(rhoecef);
    % ------------- convert to sez for calculations ---------------
    [tempvec]= rot3( rhoecef, lon          );
    [rhosez ]= rot2( tempvec, halfpi-latgd );
    %    [tempvec]= rot3( drhoecef, lon         );
    %    [drhosez]= rot2( tempvec,  halfpi-latgd);
    if rhosez(3) < 0.0
        rhosez
    end
    % ------------- calculate azimuth and elevation ---------------
    temp= sqrt( rhosez(1)*rhosez(1) + rhosez(2)*rhosez(2) );
    if ( ( temp < small ) )           % directly over the north pole
        el= sign(rhosez(3))*halfpi;   % +- 90 deg
    else
        magrhosez = mag(rhosez);
        el= asin( rhosez(3) / magrhosez );
    end

    if ( temp < small )
        az = atan2( drhosez(2), -drhosez(1) );
    else
        az= atan2( rhosez(2)/temp, -rhosez(1)/temp );
    end


