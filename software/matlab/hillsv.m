% ------------------------------------------------------------------------------
%
%                           function hillsv
%
%  this function calculates initial velocity for hills equations.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r           - initial position vector of int            m
%    alt         - altitude of tgt satellite                 km
%    dts         - desired time                              s
%
%  outputs       :
%    v           - initial velocity vector of int            m / s
%
%  locals        :
%    numkm       -
%    denom       -
%    nt          - angular velocity times time               rad
%    omega       -
%    sinnt       - sine of nt
%    cosnt       - cosine of nt
%    radius      - magnitude of range vector                 km
%
%  references    :
%    vallado       2022, 410, eq 6-60, ex 6-15
%
% [v] = hillsv( r, alt,dts );
% ------------------------------------------------------------------------------

function [v] = hillsv( r, alt,dts )

    % --------------------  implementation   ----------------------
    constastro;
    constmath;

    radius= re + alt;
    omega = sqrt( mu / (radius*radius*radius) );
    nt    = omega*dts;
    cosnt = cos( nt );
    sinnt = sin( nt );

    % --------------- determine initial velocity ------------------
    numkm= ( (6.0*r(1)*(nt-sinnt)-r(2))*omega*sinnt ...
        - 2.0*omega*r(1)*(4.0-3.0*cosnt)*(1.0-cosnt) );
    denom= (4.0*sinnt-3.0*nt)*sinnt + 4.0*( 1.0-cosnt ) ...
        *( 1.0-cosnt );

    if ( abs( denom ) > 0.000001 )
        v(2)= numkm / denom;
    else
        v(2)= 0.0;
    end
    if ( abs( sinnt ) > 0.000001 )
        v(1)= -( omega*r(1)*(4.0-3.0*cosnt) ...
            +2.0*(1.0-cosnt)*v(2) ) / ( sinnt );
    else
        v(1)= 0.0;
    end
    v(3)= -r(3)*omega*cot(nt);

end