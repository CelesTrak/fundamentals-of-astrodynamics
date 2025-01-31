% ----------------------------------------------------------------------------
%
%                           function rv2flt
%
%  this function transforms a position and velocity vector into the flight
%    elements - latgc, lon, fpa, az, position and velocity magnitude.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r           - eci position vector km
%    v           - eci velocity vector km/s
%    ttt         - julian centuries of tt                          centuries
%    jdut1       - julian date of ut1                              days from 4713 bc
%    lod         - excess length of day                            sec
%    xp          - polar motion coefficient                        arc sec
%    yp          - polar motion coefficient                        arc sec
%    terms       - number of terms for ast calculation             0,2
%    ddpsi, ddeps - corrections for fk5 to gcrf                    rad
%    iau80arr     - iau80 coefficients of eop
%
%  outputs       :
%    magr        - eci position vector magnitude                   km
%    magv        - eci velocity vector magnitude                   km/sec
%    latgc       - geocentric lat of satellite, not nadir point           -pi/2 to pi/2 rad
%    lon         - longitude                                       rad
%    fpa         - sat flight path angle                           rad
%    az          - sat flight path az                              rad
%
%  locals        :
%    fpav        - sat flight path anglefrom vert                  rad
%
%  references    :
%    vallado       2022, 111
%
% [lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt ...
%        ( reci, veci, iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps );
% ----------------------------------------------------------------------------

function [lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt ...
        ( reci, veci, iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps )
    twopi = 2.0*pi;
    small = 0.00000001;

    magr = mag(reci);
    magv = mag(veci);

    % -------- convert r to ecef for lat/lon calculation
    aeci = [0.0; 0.0; 0.0];
    [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1,...
        lod, xp, yp, 2, ddpsi, ddeps );

    % ----------------- find longitude value  ----------------- uses ecef
    temp = sqrt( recef(1)*recef(1) + recef(2)*recef(2) );
    if ( temp < small )
        lon= atan2( vecef(2), vecef(1) );
    else
        lon= atan2( recef(2), recef(1) );
    end

    %latgc = atan2( recef(3) , sqrt(recef(1)^2 + recef(2)^2) )
    latgc = asin( recef(3) / magr );

    % ------------- calculate rtasc and decl ------------------ uses eci
    temp= sqrt( reci(1)*reci(1) + reci(2)*reci(2) );
    if ( temp < small )
        rtasc= atan2( veci(2) , veci(1) );
    else
        rtasc= atan2( reci(2) , reci(1) );
    end
    %decl = atan2( reci(3) , sqrt(reci(1)^2 + reci(2)^2) )
    decl = asin( reci(3)/magr );

    h = cross(reci,veci);
    hmag = mag(h);
    rdotv= dot(reci,veci);
    fpav= atan2(hmag,rdotv);
    fpa = pi*0.5 - fpav;

    hcrossr = cross(h,reci);

    az = atan2( reci(1)*hcrossr(2) - reci(2)*hcrossr(1), hcrossr(3)*magr );

end
