%  ------------------------------------------------------------------------------
%
%                           procedure rv_elatlon
%
%  this procedure converts ecliptic latitude and longitude with position and
%    velocity vectors. uses velocity vector to find the solution of Math.Singular
%    cases.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    rijk        - ijk position vector                      km
%    vijk        - ijk velocity vector                      km / s
%    direction   - which set of vars to output              efrom  eto
%
%  outputs       :
%    rr          - radius of the sat                        km
%    ecllat      - ecliptic latitude                        -PI/2 to PI/2 rad
%    ecllon      - ecliptic longitude                       -2PI to 2PI rad
%    drr         - radius of the sat rate                   km/s
%    decllat     - ecliptic latitude rate                   -PI/2 to PI/2 rad
%    eecllon     - ecliptic longitude rate                  -2PI to 2PI rad
%
%  locals        :
%    obliquity   - obliquity of the ecliptic                rad
%    temp        -
%    temp1       -
%    re          - position vec in eclitpic frame
%    ve          - velocity vec in ecliptic frame
%
%  coupling      :
%    mag         - magnitude of a vector
%    rot1        - rotation about 1st axis
%    dot         - dot product
%    arcsin      - arc Math.Sine function
%    Math.Atan2       - arc tangent function that resolves quadrant ambiguites
%
%  references    :
%    vallado       2022 268, eq 4-15
%
% [rr, ecllon, ecllat, drr, decllon, decllat] = rv2elatlon (rijk,vijk);
% ----------------------------------------------------------------------------

function [rr, ecllon, ecllat, drr, decllon, decllat] = rv2elatlon (rijk,vijk)
    % --------------------  implementation   ----------------------
    constmath;
    obliquity= 0.40909280;   %23.439291 /rad

    [r] = rot1 ( rijk, obliquity );
    [v] = rot1 ( vijk, obliquity );

    % ------------- calculate angles and rates ----------------
    rr= mag(r);
    temp= sqrt( r(1)*r(1) + r(2)*r(2) );
    if ( temp < small )
        temp1= sqrt( v(1)*v(1) + v(2)*v(2) );
        if ( abs(temp1) > small )
            ecllon= atan2( v(2) , v(1) );
        else
            ecllon= 0.0;
        end
    else
        ecllon= atan2( r(2) , r(1) );
    end
    ecllat= asin( r(3)/rr );

    temp1= -r(2)*r(2) - r(1)*r(1);  % different now
    drr= dot(r,v)/rr;
    if ( abs( temp1 ) > small )
        decllon= ( v(1)*r(2) - v(2)*r(1) ) / temp1;
    else
        decllon= 0.0;
    end
    if ( abs( temp ) > small )
        decllat= ( v(3) - drr*sin( ecllat ) ) / temp;
    else
        decllat= 0.0;
    end

end