% ------------------------------------------------------------------------------
%
%                            function radec2rv
%
%   this function converts the right ascension and declination values with
%     position and velocity vectors of a satellite. uses velocity vector to
%     find the solution of singular cases.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%     r           - position vector eci                       km
%     v           - velocity vector eci                       km/s
%     direct      -  direction to convert                     eFrom  eTo
%
%   outputs       :
%     rr          - radius of the satellite                      km
%     rtasc       - right ascension                              rad
%     decl        - declination                                  rad
%     drr         - radius of the satellite rate                 km/s
%     drtasc      - right ascension rate                         rad/s
%     ddecl       - declination rate                             rad/s
%
%   locals        :
%     temp        - temporary position vector
%     temp1       - temporary variable
%
%   coupling      :
%     none
%
%   references    :
%    vallado       2022, 256, alg 25
%
% [r, v] = radec2rv( rr, rtasc, decl, drr, drtasc, ddecl );
% ------------------------------------------------------------------------------

function [r, v] = radec2rv( rr, rtasc, decl, drr, drtasc, ddecl )
    % -------------------------  implementation   -----------------
    small        = 0.00000001;

    r(1)= rr*cos(decl)*cos(rtasc);
    r(2)= rr*cos(decl)*sin(rtasc);
    r(3)= rr*sin(decl);
    r = r';

    v(1)= drr*cos(decl)*cos(rtasc) - rr*sin(decl)*cos(rtasc)*ddecl ...
        - rr*cos(decl)*sin(rtasc)*drtasc;
    v(2)= drr*cos(decl)*sin(rtasc) - rr*sin(decl)*sin(rtasc)*ddecl ...
        + rr*cos(decl)*cos(rtasc)*drtasc;
    v(3)= drr*sin(decl) + rr*cos(decl)*ddecl;
    v = v';
end