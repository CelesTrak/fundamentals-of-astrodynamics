% -----------------------------------------------------------------------------
%
%                           function lon2nu
%
%  this function finds the true anomaly from coes and longitude.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    jdut1       - julian date in ut1                       days from 4713 bc
%    lon         - longitude                                0 to 2pi rad
%    incl        - inclination                              0 to 2pi rad
%    raan        - right ascenion of the node               0 to 2pi rad
%    argp        - argument of perigee                      0 to 2pi rad
%
%  outputs       :
%    nu          - true anomaly                             0 to 2pi rad
%
%  locals        :
%    temp        - temporary variable for doubles           rad
%    tut1        - julian centuries from the
%                  jan 1, 2000 12 h epoch (ut1)
%
%  coupling      :
%    none
%
%  references    :
%    vallado       2022, 112, eq 2-103
%
%  nu = lon2nu( jdut1, lon, incl, raan, argp);
% -----------------------------------------------------------------------------

function nu = lon2nu( jdut1, lon, incl, raan, argp)
    rad = 180/pi;
    twopi = 2.0 * pi;

    %    fprintf(' jd %16.8f lon %11.5f  incl %11.5f raan %11.5f argp %11.5f \n',jdut1, lon*rad, incl*rad, raan*rad, argp*rad );
    %    need to use their GMST calculation
    ed = jdut1 +0.0 - 2451544.5;  % elapsed days from 1 jan 2000 0 hrs
    gmst = 99.96779469 + 360.9856473662860 * ed + 0.29079e-12 * ed * ed;  % deg
    deg2rad    = pi/180.0;
    gmst = rem( gmst*deg2rad,2.0*pi );

    % ------------------------ check quadrants --------------------
    if ( gmst < 0.0 )
        gmst = gmst + 2.0*pi;
    end

    lambdau = gmst + lon - raan;

    % make sure lambdau is 0 to 360 deg
    if lambdau < 0.0
        lambdau = lambdau + 2.0*pi;
    end
    if lambdau > twopi
        lambdau = lambdau - 2.0*pi;
    end

    arglat = atan(tan(lambdau) / cos(incl));
    % find nu
    if (lambdau >= 0.5*pi) && (lambdau <  1.5*pi)
        arglat = arglat + pi;
    end

    temp = arglat - argp;

    nu = temp;
    % fprintf(' %11.5f %11.5f %11.5f %11.5f  %16.10f ',lambdau*rad, argp*rad, lon*rad, gmst*rad, temp*rad);
    fprintf(' lu %11.5f argp %11.5f lon %11.5f gmst %11.5f arglat %11.5f nu %16.10f ',lambdau*rad, argp*rad, lon*rad, gmst*rad, arglat*rad, temp*rad);

    %     fprintf(' nu = %11.5f deg \n',nu*rad);
end