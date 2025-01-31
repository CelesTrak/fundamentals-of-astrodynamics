% -----------------------------------------------------------------------------
%
%                           function rv2pqw
%
%  this function finds the pqw vectors given the geocentric equatorial
%  position and velocity vectors.  mu is needed if km and m are
%    both used with the same routine
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r           - ijk position vector            km or m
%    v           - ijk velocity vector            km/s or m/s
%    mu          - gravitational parameter        km3/s2 or m3/s2
%
%  outputs       :
%    rpqw        - pqw position vector            km
%    vpqw        - pqw velocity vector            km / s
%
%  locals        :
%    hbar        - angular momentum h vector      km2 / s
%    ebar        - eccentricity     e vector
%    nbar        - line of nodes    n vector
%    c1          - v**2 - u/r
%    rdotv       - r dot v
%    hk          - hk unit vector
%    sme         - specfic mechanical energy      km2 / s2
%    i           - index
%    p           - semilatus rectum               km
%    a           - semimajor axis                 km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to Math.PI rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
%    truelon     - true longitude            (ce) 0.0  to 2pi rad
%    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
%    temp        - temporary variable
%    typeorbit   - type of orbit                  ee, ei, ce, ci
%
%  coupling      :
%    mag         - magnitude of a vector
%    MathTimeLibr.cross       - MathTimeLibr.cross product of two vectors
%    angle       - find the angle between two vectors
%
%  references    :
%    vallado       2022, 126, alg 9, ex 2-5
%
% [rpqw, vpqw] = rv2pqw(r, v);
% -----------------------------------------------------------------------------

function [rpqw, vpqw] = rv2pqw(r, v)
    constastro;
    small = 0.00000001;
    undefined = 999999.1;

    [p,a,ecc,incl,raan,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r,v);

    if (nu < undefined)
        % ----------  form pqw position and velocity vectors ----------
        sin_nu = sin(nu);
        cos_nu = cos(nu);
        temp = p / (1.0 + ecc * cos_nu);
        rpqw(1) = temp * cos_nu;
        rpqw(2) = temp * sin_nu;
        rpqw(3) = 0.0;
        if (abs(p) < 0.00000001)
            p = 0.00000001;
        end

        vpqw(1) = -sin_nu * sqrt(mu / p);
        vpqw(2) = (ecc + cos_nu) * sqrt(mu / p);
        vpqw(3) = 0.0;
    else
        rpqw(1) = undefined;
        rpqw(2) = undefined;
        rpqw(3) = undefined;
        vpqw(1) = undefined;
        vpqw(2) = undefined;
        vpqw(3) = undefined;
    end
end  % rv2pqw