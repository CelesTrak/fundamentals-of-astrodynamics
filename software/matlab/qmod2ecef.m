
% ----------------------------------------------------------------------------
%
%                           function qmod2ecef
%
%  this function trsnforms a vector from the mean equator mean equniox frame
%    (j2000), to an earth fixed (ITRF) frame.  the results take into account
%    the effects of precession, nutation, sidereal time, and polar motion.
%
%  author        : david vallado                  719-573-2600   27 jun 2002
%
%  revisions
%    vallado     - add terms for ast calculation                 30 sep 2002
%    vallado     - consolidate with iau 2000                     14 feb 2005
%
%  inputs          description                    range / units
%    reci        - position vector eci            km
%    veci        - velocity vector eci            km/s
%    aeci        - acceleration vector eci        km/s2
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    eqeterms    - terms for ast calculation      0,2
%    ddpsi       - delta psi correction to gcrf   rad
%    ddeps       - delta eps correction to gcrf   rad
%
%  outputs       :
%    recef       - position vector earth fixed    km
%    vecef       - velocity vector earth fixed    km/s
%    aecef       - acceleration vector earth fixedkm/s2
%
%  locals        :
%    deltapsi    - nutation angle                 rad
%    trueeps     - true obliquity of the ecliptic rad
%    meaneps     - mean obliquity of the ecliptic rad
%    omega       -                                rad
%    prec        - matrix for mod - eci
%    nut         - matrix for tod - mod
%    st          - matrix for pef - tod
%    stdot       - matrix for pef - tod rate
%    pm          - matrix for ecef - pef
%
%  coupling      :
%   precess      - rotation for precession        eci - mod
%   nutation     - rotation for nutation          mod - tod
%   sidereal     - rotation for sidereal time     tod - pef
%   polarm       - rotation for polar motion      pef - ecef
%
%  references    :
%    vallado       2004, 219-228
%
% [recef,vecef,aecef] = qmod2ecef  ( rqmod, vqmod, iau80arr, ttt, jdutc )
% ----------------------------------------------------------------------------

function [recef,vecef,aecef] = qmod2ecef  ( rqmod, vqmod, iau80arr, ttt, jdutc )

    [deltapsi, trueeps, meaneps, omega, nut] = nutationqmod(iau80arr, ttt);

    %gmstx= gstime( jdutc );

    ed= jdutc - 2451544.5;  % elapsed days from 1 jan 2000 0 hr
    gmst = 99.96779469 + 360.9856473662860 * ed + 0.29079e-12 * ed * ed;  % deg
    deg2rad    = pi/180.0;
    gmst = rem( gmst*deg2rad,2.0*pi );

    % ------------------------ check quadrants --------------------
    if ( gmst < 0.0 )
        gmst = gmst + 2.0*pi;
    end

    % do ast and see if it affects
    if (jdutc > 2450449.5 )
        ast= gmst + deltapsi* cos(meaneps) ...
            + 0.00264*pi /(3600*180)*sin(omega) ...
            + 0.000063*pi /(3600*180)*sin(2.0 *omega);
    else
        ast= gmst + deltapsi* cos(meaneps);
    end
    ast = rem (ast,2*pi);
    %gmst = ast;


    thetasa    = 7.29211514670698e-05; %  * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    st(1,1) =  cos(gmst);
    st(1,2) = -sin(gmst);
    st(1,3) =  0.0;
    st(2,1) =  sin(gmst);
    st(2,2) =  cos(gmst);
    st(2,3) =  0.0;
    st(3,1) =  0.0;
    st(3,2) =  0.0;
    st(3,3) =  1.0;

    recef  = st'*nut'*reci;

    vecef  = st'*nut'*veci - cross( omegaearth,recef );

    temp  = cross(omegaearth,recef);

    aecef = st'*nut'*aeci - cross(omegaearth,temp) - 2.0*cross(omegaearth,vecef);

