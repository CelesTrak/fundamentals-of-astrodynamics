% ------------------------------------------------------------------------------
%
%                           procedure rv2razel
%
%  this procedure converts range, azimuth, and elevation and their rates with
%    the geocentric equatorial (ecef) position and velocity vectors.  notice the
%    value of small as it can affect rate term calculations. uses velocity
%    vector to find the solution of singular cases.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    recef       - ecef position vector                       km
%    vecef       - ecef velocity vector                       km/s
%    latgd       - geodetic latitude                          -pi/2 to pi/2 rad
%    lon         - geodetic longitude                         -2pi to pi rad
%    direct      -  direction to convert                      eFrom  eTo
%
%  outputs       :
%    rho         - satellite range from site                  km
%    az          - azimuth                                    0.0 to 2pi rad
%    el          - elevation                                  -pi/2 to pi/2 rad
%    drho        - range rate                                 km/s
%    daz         - azimuth rate                               rad/s
%    del         - elevation rate                             rad/s
%
%  locals        :
%    rsecef      - ecef site position vector                  km
%    rhovecef    - ecef range vector from site                km
%    drhovecef   - ecef velocity vector from site             km/s
%    rhosez      - sez range vector from site                 km
%    drhosez     - sez velocity vector from site              km
%    tempvec     - temporary vector
%    temp        - temporary extended value
%    temp1       - temporary extended value
%    i           - index
%
%  coupling      :
%    mag         - magnitude of a vector
%    addvec      - add two vectors
%    rot3        - rotation about the 3rd axis
%    rot2        - rotation about the 2nd axis
%    atan2       - arc tangent function which also resloves quadrants
%    dot         - dot product of two vectors
%    rvsez_razel - find r and v from site in topocentric horizon (sez) system
%    arcsin      - arc sine function
%    sign        - returns the sign of a variable
%
%  references    :
%    vallado       2022, 262, alg 27
%
% [rho, az, el, drho, daz, del] = rv2razel(recef, vecef, latgd, lon, alt );
% ------------------------------------------------------------------------------

function [rho, az, el, drho, daz, del] = rv2razel(recef, vecef, latgd, lon, alt )
    halfpi = pi*0.5;
    small  = 0.00000001;

    % --------------------- implementation ------------------------
    % ----------------- get site vector in ecef -------------------
    [rsecef, vsecef] = site (latgd, lon, alt );
    %fprintf(1,'rsecef    %14.7f %14.7f %14.7f \n',rsecef );

    % -------------------- convert eci to ecef --------------------
    % aeci = [0;0;0];
    % [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    %fprintf(1,'sat recef    %14.7f %14.7f %14.7f \n',recef );
    % simplified - just use sidereal time rotation
    % thetasa= earthrot * (1.0  - 0.0/86400.0 );
    % omegaearth = [0; 0; thetasa;];
    % [deltapsi,trueeps,meaneps,omega,nut] = nutation(ttt,ddpsi,ddeps);
    % [st,stdot] = sidereal(jdut1,deltapsi,meaneps,omega,0,0 );
    %  recef  = st'*reci;
    %  vecef  = st'*veci - cross( omegaearth,recef );


    % ------- find ecef range vector from site to satellite -------
    rhoecef  = recef - rsecef;
    drhoecef = vecef;
    rho      = mag(rhoecef);

    % ------------- convert to sez for calculations ---------------
    [tempvec]= rot3( rhoecef, lon          );
    [rhosez ]= rot2( tempvec, halfpi-latgd );

    [tempvec]= rot3( drhoecef, lon         );
    [drhosez]= rot2( tempvec,  halfpi-latgd);

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

    % ------ calculate range, azimuth and elevation rates ---------
    drho= dot(rhosez,drhosez)/rho;
    if ( abs( temp*temp ) > small )
        daz= ( drhosez(1)*rhosez(2) - drhosez(2)*rhosez(1) ) / ( temp*temp );
    else
        daz= 0.0;
    end

    if ( abs( temp ) > small )
        del= ( drhosez(3) - drho*sin( el ) ) / temp;
    else
        del= 0.0;
    end

end