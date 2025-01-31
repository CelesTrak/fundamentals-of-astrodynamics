%  ---------------------------------------------------------------------------
%
%                           function gd2gc
%
%  this function converts from geodetic to geocentric latitude for positions
%    on the surface of the earth.  notice that (1-f) squared = 1-esqrd.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    latgd       - geodetic latitude                        -PI to PI rad
%
%  outputs       :
%    latgc       - geocentric lat of satellite, not nadir point  -pi/2 to pi/2 rad
%
%  locals        :
%    none.
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2022, 146, eq 3-14
%
% [latgc] = gd2gc ( latgd );
% ------------------------------------------------------------------------------

function [latgc] = gd2gc ( latgd )

    eesqrd = 0.006694385000;     % eccentricity of earth sqrd

    % -------------------------  implementation   -----------------
    latgc= atan( (1.0  - eesqrd)*tan(latgd) );

end