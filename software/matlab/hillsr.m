% ------------------------------------------------------------------------------
%
%                           function hillsr
%
%  this function calculates various position information for hills equations.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r           - init rel position of int                 m or km
%    v           - init rel velocity of int                 m or km/s
%    alt         - altitude of tgt satellite                km
%    dts         - desired time                             s
%
%  outputs       :
%    rinit       - final rel position of int                m or km
%    vinit       - final rel velocity of int                m or km/s
%
%  locals        :
%    nt          - angular velocity times time              rad
%    omega       -
%    sinnt       - sine of nt
%    cosnt       - cosine of nt
%    radius      - magnitude of range vector                km
%
%  references    :
%    vallado       2022, 401, alg 48, ex 6-14
%
%  [rint, vint] = hillsr( r,v, alt, dts );
% ------------------------------------------------------------------------------

function [rint, vint] = hillsr( r, v, alt, dts )

    % --------------------  implementation   ----------------------
    constastro;

    radius= re + alt; % in km
    omega = sqrt( mu / (radius*radius*radius) ); % rad/s
    nt    = omega * dts;
    cosnt = cos( nt );
    sinnt = sin( nt );

    % --------------- determine new positions  --------------------
    rint(1)= ( v(1)/omega ) * sinnt - ...
        ( (2.0*v(2)/omega) + 3.0*r(1) ) * cosnt + ...
        ( (2.0*v(2)/omega) + 4.0*r(1) );
    rint(2)= ( 2.0*v(1)/omega ) * cosnt + ...
        ( (4.0*v(2)/omega) + 6.0*r(1) ) * sinnt + ...
        ( r(2) - (2.0*v(1)/omega) ) - ...
        ( 3.0*v(2) + 6.0*omega*r(1) )*dts;
    rint(3)= r(3)*cosnt + (v(3)/omega)*sinnt;

    % --------------- determine new velocities  -------------------
    vint(1)= v(1)*cosnt + (2.0*v(2)+3.0*omega*r(1))*sinnt;
    vint(2)= -2.0*v(1)*sinnt + (4.0*v(2) ...
        +6.0*omega*r(1))*cosnt - (3.0*v(2)+6.0*omega*r(1));
    vint(3)= -r(3)*omega*sinnt + v(3)*cosnt;

end
