% ----------------------------------------------------------------------------
%
%                           function flt2rv.m
%
%  this function transforms  the flight elements - latgc, lon, fpav, az,
%    position and velocity magnitude into an eci position and velocity vector.
%
%  author        : david vallado                  719-573-2600   17 jun 2002
%
%  revisions
%    vallado     - fix extra terms in rtasc calc                  8 oct 2002
%
%  inputs          description                    range / units
%    rmag        - eci position vector magnitude  km
%    vmag        - eci velocity vector magnitude  km/sec
%    latgc       - geocentric latitude            rad
%    lon         - longitude                      rad
%    fpa         - sat flight path angle          rad
%    az          - sat flight path az             rad
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    terms       - number of terms for ast calculation 0,2
%    ddpsi,ddeps - corrections for fk5 to gcrf    rad
%
%  outputs       :
%    r           - eci position vector            km
%    v           - eci velocity vector            km/s
%
%  locals        :
%    fpav        - sat flight path anglefrom vert rad
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado       2013, xx
%    escobal            397
%    chobotov            67
%
% [reci, veci] = flt2rv ( rmag,vmag,latgc,lon,fpa,az,iau80arr,ttt,jdut1,lod,xp,yp,ddpsi,ddeps );
% ----------------------------------------------------------------------------

function [reci, veci] = flt2rv ( rmag,vmag,latgc,lon,fpa,az,iau80arr, ttt,jdut1,lod,xp,yp,ddpsi,ddeps )

    twopi = 2.0*pi;
    small        = 0.00000001;

    % -------- form position vector
    recef(1) = rmag*cos(latgc)*cos(lon);
    recef(2) = rmag*cos(latgc)*sin(lon);
    recef(3) = rmag*sin(latgc);
    recef=recef';

    % -------- convert r to eci
    vecef = [0.0; 0.0; 0.0];  % this is a dummy for now
    aecef = [0.0; 0.0; 0.0];
    [reci, veci, aeci] = ecef2eci(recef, vecef, aecef, iau80arr, ttt, ...
        jdut1, lod, xp, yp, 2, ddpsi, ddeps );

    % ------------- calculate rtasc and decl ------------------
    temp= sqrt( reci(1)*reci(1) + reci(2)*reci(2) );

    if ( temp < small )
        % v needs to be defined herexxxxxxxxx
        rtasc= atan2( veci(2) , veci(1) );
    else
        rtasc= atan2( reci(2) , reci(1) );
    end
    decl= asin( reci(3)/rmag );

    % -------- form velocity vector
    fpav = pi*0.5 - fpa;
    veci(1)= vmag*( -cos(rtasc)*sin(decl)*(cos(az)*cos(fpav) - ...
        sin(rtasc)*sin(az)*cos(fpav)) + cos(rtasc)*sin(decl)*sin(fpav) );
    veci(2)= vmag*( -sin(rtasc)*sin(decl)*(cos(az)*cos(fpav) + ...
        cos(rtasc)*sin(az)*cos(fpav)) + sin(rtasc)*cos(decl)*sin(fpav) );
    veci(3)= vmag*(  sin(decl)*sin(fpav) + cos(decl)*cos(az)*cos(fpav) );

end