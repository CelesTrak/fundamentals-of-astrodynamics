
% ------------------------------------------------------------------------------
%
%                           procedure razel2rv
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
% [recef, vecef] = razel2rv(latgd, lon, alt, rho, az, el, drho, daz, del);
% ------------------------------------------------------------------------------

function [recef, vecef] = razel2rv(latgd, lon, alt, rho, az, el, drho, daz, del)

        % -------------------------  implementation   -----------------
        constmath;

        % -----------  find sez range and velocity vectors ------------
        [rhosez,drhosez] = raz2rvs( rho,az,el,drho,daz,del );

        % -----------  perform sez to ijk (ecef) transformation -------
        [tempvec] = rot2( rhosez , latgd-halfpi );
        [rhoecef] = rot3( tempvec,-lon          );
        rhoecef = rhoecef';
        
        [tempvec] = rot2( drhosez, latgd-halfpi );
        [drhoecef]= rot3( tempvec,-lon          );
        drhoecef = drhoecef';

        % ----------  find ecef range and velocity vectors -------------
        [rs, vs] = site ( latgd, lon, alt );
        recef = rhoecef + rs;
        vecef = drhoecef;

        % % -------- convert ecef to eci
        % recef = recef;
        % vecef = vecef;
        % aceef     = [0;0;0];
        % [reci,veci,aeci] = ecef2eci(recef,vecef,acef,ttt,iau80arr, jdut1,lod,xp,yp,terms,ddpsi,ddeps );

end