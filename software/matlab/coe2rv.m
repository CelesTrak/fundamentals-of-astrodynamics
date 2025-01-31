% ------------------------------------------------------------------------------
%
%                           function coe2rv
%
%  this function finds the position and velocity vectors in geocentric
%    equatorial (ijk) system given the classical orbit elements. the additional
%    orbital elements provide calculations for perfectly circular and equatorial orbits.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    p           - semilatus rectum                          km
%    ecc         - eccentricity
%    incl        - inclination                               0.0 to pi rad
%    raan        - rtasc of ascending node                   0.0 to 2pi rad
%    argp        - argument of perigee                       0.0 to 2pi rad
%    nu          - true anomaly                              0.0 to 2pi rad
%    arglat      - argument of latitude                      (ci) 0.0 to 2pi rad
%    lamtrue     - true longitude                            (ce) 0.0 to 2pi rad
%    lonper      - longitude of periapsis                    (ee) 0.0 to 2pi rad
%
%  outputs       :
%    r           - ijk position vector                        km
%    v           - ijk velocity vector                        km / s
%
%  locals        :
%    temp        - temporary real*8 value
%    rpqw        - pqw position vector                        km
%    vpqw        - pqw velocity vector                        km / s
%    sinnu       - sine of nu
%    cosnu       - cosine of nu
%    tempvec     - pqw velocity vector
%
%  coupling      :
%    rot3        - rotation about the 3rd axis
%    rot1        - rotation about the 1st axis
%
%  references    :
%    vallado       2022, 120, alg 10, ex 2-5
%
% [r,v] = coe2rv ( p,ecc,incl,omega,argp,nu,arglat,truelon,lonper );
% ------------------------------------------------------------------------------

function [r,v] = coe2rv ( p,ecc,incl,raan,argp,nu,arglat,truelon,lonper )

    % -------------------------  implementation   -----------------
    constmath;
    constastro;

    % -------------------------------------------------------------
    %       determine what type of orbit is involved and set up the
    %       set up angles for the special cases.
    % -------------------------------------------------------------
    if ( ecc < small )
        % ----------------  circular equatorial  ------------------
        if (incl<small) || ( abs(incl-pi)< small )
            argp = 0.0;
            raan= 0.0;
            nu   = truelon;
        else
            % --------------  circular inclined  ------------------
            argp= 0.0;
            nu  = arglat;
        end
    else
        % ---------------  elliptical equatorial  -----------------
        if ( ( incl<small) || (abs(incl-pi)<small) )
            argp = lonper;
            raan= 0.0;
        end
    end

    % ----------  form pqw position and velocity vectors ----------
    cosnu = cos(nu);
    sinnu = sin(nu);
    temp  = p / (1.0  + ecc*cosnu);
    rpqw(1) = temp*cosnu;
    rpqw(2) = temp*sinnu;
    rpqw(3) =     0.0;
    if ( abs(p) < 0.0001)
        p = 0.0001;
    end
    vpqw(1) =    -sinnu*sqrt(mu)  / sqrt(p);
    vpqw(2) =  (ecc + cosnu)*sqrt(mu) / sqrt(p);
    vpqw(3) =      0.0;

    % ----------------  perform transformation to ijk  ------------
    [tempvec] = rot3( rpqw   , -argp );
    [tempvec] = rot1( tempvec, -incl );
    [r] = rot3( tempvec, -raan );

    [tempvec] =rot3( vpqw   , -argp );
    [tempvec] =rot1( tempvec, -incl );
    [v] = rot3( tempvec, -raan );

    r = r';
    v = v';

end