%  ---------------------------------------------------------------------------
%
%                           procedure site
%
%  this function finds the position and velocity vectors for a site.  the
%    answer is returned in the geocentric equatorial (ecef) coordinate system.
%    note that the velocity is zero because the coordinate system is fixed to
%    the earth.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    latgd       - geodetic latitude                        -PI/2 to PI/2 rad
%    lon         - longitude of site                        -2pi to 2pi rad
%    alt         - altitude                                 km
%
%  outputs       :
%    rsecef      - ecef site position vector                km
%    vsecef      - ecef site velocity vector                km/s
%
%  locals        :
%    sinlat      - variable containing sin(lat)             rad
%    temp        - temporary real value
%    rdel        - rdel component of site vector            km
%    rk          - rk component of site vector              km
%    cearth      -
%
%  references    :
%    vallado       2022, 436, alg 51, ex 7-1
%
% [rsecef, vsecef] = site ( latgd,lon,alt );
% -----------------------------------------------------------------------------

function [rsecef, vsecef] = site ( latgd,lon,alt )
    constastro;

    % -------------------------  implementation   -----------------
    sinlat      = sin( latgd );

    % ------  find rdel and rk components of site vector  ---------
    cearth= re / sqrt( 1.0 - ( eccearthsqrd*sinlat*sinlat ) );
    rdel  = ( cearth + alt )*cos( latgd );
    rk    = ( (1.0-eccearthsqrd)*cearth + alt )*sinlat;

    % ---------------  find site position vector  -----------------
    rsecef(1) = rdel * cos( lon );
    rsecef(2) = rdel * sin( lon );
    rsecef(3) = rk;
    rsecef = rsecef';

    % ---------------  find site velocity vector  -----------------
    vsecef = [0.0; 0.0; 0.0];

end