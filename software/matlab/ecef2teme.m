% ----------------------------------------------------------------------------
%
%                           function ecef2teme
%
%  this function trsnforms a vector from the earth fixed (ITRF) frame to the
%    true equator mean equniox frame (teme). the results take into account
%    the effects of sidereal time, and polar motion.
%
%  author        : david vallado                  719-573-2600   30 oct 2017
%
%  revisions
%
%  inputs          description                    range / units
%    recef       - position vector earth fixed    km
%    vecef       - velocity vector earth fixed    km/s
%    aecef       - acceleration vector earth fixedkm/s2
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    eqeterms    - use extra two terms (kinematic) after 1997  0, 2
%
%  outputs       :
%    rteme       - position vector teme           km
%    vteme       - velocity vector teme           km/s
%    ateme       - acceleration vector teme       km/s2
%
%  locals        :
%    st          - matrix for pef - tod
%    pm          - matrix for ecef - pef
%
%  coupling      :
%   gstime       - greenwich mean sidereal time   rad
%   polarm       - rotation for polar motion      pef - ecef
%
%  references    :
%    vallado       2013, 231-233
%
% [rteme, vteme, ateme] = ecef2teme( recef, vecef, aecef, ttt, jdut1, lod, xp, yp, eqeterms );
% ----------------------------------------------------------------------------

function [rteme, vteme, ateme] = ecef2teme( recef, vecef, aecef, ttt, jdut1, lod, xp, yp, eqeterms )
    constastro;
    deg2rad = pi/180.0;

    % ------------------------ find gmst --------------------------
    gmst= gstime( jdut1 );

    % find omeage from nutation theory
    omega=  125.04452222  + (   -6962890.5390 *ttt + ...
        7.455 *ttt*ttt + 0.008 *ttt*ttt*ttt )  / 3600.0;
    omega= rem( omega, 360.0  ) * deg2rad;

    % teme does not include the geometric terms here
    % after 1997, kinematic terms apply
    if (jdut1 > 2450449.5 ) && (eqeterms > 0)
        gmstg = gmst ...
            + 0.00264*pi /(3600*180)*sin(omega) ...
            + 0.000063*pi /(3600*180)*sin(2.0 *omega);
    else
        gmstg = gmst;
    end

    gmstg = rem (gmstg, 2.0*pi);


    st(1,1) =  cos(gmstg);
    st(1,2) = -sin(gmstg);
    st(1,3) =  0.0;
    st(2,1) =  sin(gmstg);
    st(2,2) =  cos(gmstg);
    st(2,3) =  0.0;
    st(3,1) =  0.0;
    st(3,2) =  0.0;
    st(3,3) =  1.0;

    [pm] = polarm(xp, yp, ttt, '80');

    rpef = pm*recef;
    rteme = st*rpef;

    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    vpef = pm*vecef;
    vteme = st*(vpef + cross(omegaearth,rpef));

    temp = cross(omegaearth,rpef);
    ateme = st*( pm*aecef + cross(omegaearth,temp) ...
        + 2.0*cross(omegaearth,vpef) );


