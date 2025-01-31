% -----------------------------------------------------------------------------
%
%                           function rv2coe
%
%  this function finds the classical orbital elements given the geocentric
%    equatorial position and velocity vectors. mu is needed if km and m are
%    both used with the same routine
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r           - ijk position vector                          km
%    v           - ijk velocity vector                          km/s
%
%  outputs       :
%    p           - semilatus rectum                             km
%    a           - semimajor axis                               km
%    ecc         - eccentricity
%    incl        - inclination                                  0.0  to Math.PI rad
%    raan       - right ascension of ascending node             0.0  to 2pi rad
%    argp        - argument of perigee                          0.0  to 2pi rad
%    nu          - true anomaly                                 0.0  to 2pi rad
%    m           - mean anomaly                                 0.0  to 2pi rad
%    arglat      - argument of latitude      (ci)               0.0  to 2pi rad
%    truelon     - true longitude            (ce)               0.0  to 2pi rad
%    lonper      - longitude of periapsis    (ee)               0.0  to 2pi rad
%
%  locals        :
%    hbar        - angular momentum h vector                    km2 / s
%    ebar        - eccentricity     e vector
%    nbar        - line of nodes    n vector
%    c1          - v**2 - u/r
%    rdotv       - r dot v
%    hk          - hk unit vector
%    sme         - specfic mechanical energy                    km2 / s2
%    i           - index
%    e           - eccentric, parabolic,
%                  hyperbolic anomaly                           rad
%    temp        - temporary variable
%    typeorbit   - type of orbit                                ee, ei, ce, ci
%
%  coupling      :
%    mag         - magnitude of a vector
%    cross       - cross product of two vectors
%    angle       - find the angle between two vectors
%    newtonnu    - find the mean anomaly
%
%  references    :
%    vallado       2022, 115, alg 9, ex 2-5
%
% [p, a, ecc, incl, raan, argp, nu, m, arglat,truelon,lonper] = rv2coe (r, v);
% ------------------------------------------------------------------------------

function [p, a, ecc, incl, raan, argp, nu, m, arglat,truelon,lonper] = rv2coe (r, v)
    constmath;
    constastro;  % don't overwrite mu
    muin  = mu;  % this is the km version
    small = 1e-12;

    % -------------------------  implementation   -----------------
    magr = mag( r );
    magv = mag( v );
    % ------------------  find h n and e vectors   ----------------
    [hbar] = cross( r,v );
    magh = mag( hbar );
    if ( magh >= 0.0 )
        nbar(1)= -hbar(2);
        nbar(2)= hbar(1);
        nbar(3)= 0.0;
        magn = mag( nbar );
        c1 = magv*magv - muin /magr;
        rdotv= dot( r,v );
        for i= 1 : 3
            ebar(i)= (c1*r(i) - rdotv*v(i))/muin;
        end
        ecc = mag( ebar );

        % ------------  find a e and semi-latus rectum   ----------
        sme= ( magv*magv*0.5  ) - ( muin /magr );
        if ( abs( sme ) > small )
            a= -muin  / (2.0 *sme);
        else
            a= infinite;
        end
        p = magh*magh/muin;

        % -----------------  find inclination   -------------------
        hk= hbar(3)/magh;
        incl= acos( hk );

        % --------  determine type of orbit for later use  --------
        % ------ elliptical, parabolic, hyperbolic inclined -------
        typeorbit= 'ei';
        if ( ecc < small )
            % ----------------  circular equatorial ---------------
            if  (incl<small) || (abs(incl-pi)<small)
                typeorbit= 'ce';
            else
                % --------------  circular inclined ---------------
                typeorbit= 'ci';
            end
        else
            % - elliptical, parabolic, hyperbolic equatorial --
            if  (incl<small) || (abs(incl-pi)<small)
                typeorbit= 'ee';
            end
        end

        % ----------  find right ascension of ascending node ------------
        if ( magn > small )
            temp= nbar(1) / magn;
            if ( abs(temp) > 1.0  )
                temp= sign(temp);
            end
            raan= acos( temp );
            if ( nbar(2) < 0.0  )
                raan= twopi - raan;
            end
        else
            raan= undefined;
        end

        % ---------------- find argument of perigee ---------------
        if (strcmp(typeorbit, 'ei') == 1)
            argp = angl( nbar,ebar);
            if ( ebar(3) < 0.0  )
                argp= twopi - argp;
            end
        else
            argp= undefined;
        end

        % ------------  find true anomaly at epoch    -------------
        if ( typeorbit(1:1) == 'e' )
            nu =  angl( ebar,r);
            if ( rdotv < 0.0  )
                nu= twopi - nu;
            end
        else
            nu= undefined;
        end

        % ----  find argument of latitude - circular inclined -----
        % -- find in general cases too
        if (strcmp(typeorbit, 'ci') == 1) || (strcmp(typeorbit, 'ei') == 1)
            arglat = angl( nbar,r );
            if ( r(3) < 0.0  )
                arglat= twopi - arglat;
            end
            m = arglat;
        else
            arglat= undefined;
        end

        % -- find longitude of perigee - elliptical equatorial ----
        if  ( ecc>small ) && (strcmp(typeorbit, 'ee') == 1)
            temp= ebar(1)/ecc;
            if ( abs(temp) > 1.0  )
                temp= sign(temp);
            end
            lonper= acos( temp );
            if ( ebar(2) < 0.0  )
                lonper= twopi - lonper;
            end
            if ( incl > halfpi )
                lonper= twopi - lonper;
            end
        else
            lonper= undefined;
        end

        % -------- find true longitude - circular equatorial ------
        if  ( magr>small ) && (strcmp(typeorbit, 'ce') == 1)
            temp= r(1)/magr;
            if ( abs(temp) > 1.0  )
                temp= sign(temp);
            end
            truelon= acos( temp );
            if ( r(2) < 0.0  )
                truelon= twopi - truelon;
            end
            if ( incl > halfpi )
                truelon= twopi - truelon;
            end
            m = truelon;
        else
            truelon= undefined;
        end

        % ------------ find mean anomaly for all orbits -----------
        if (isreal(ecc) )
            [e,m] = newtonnu(ecc, nu);
        end

    else
        p    = undefined;
        a    = undefined;
        ecc  = undefined;
        incl = undefined;
        raan= undefined;
        argp = undefined;
        nu   = undefined;
        m    = undefined;
        arglat = undefined;
        truelon= undefined;
        lonper = undefined;
    end

end