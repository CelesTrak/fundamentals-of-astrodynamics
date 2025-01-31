% ----------------------------------------------------------------------------
%
%                           function mod_eci
%
%  this function transforms between the mean equator mean equinox (eci), and
%    the mean of date frame (mod).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    reci        - position vector eci                           km
%    veci        - velocity vector eci                           km/s
%    direct      - direction                                     eto, efrom
%    iau80arr    - iau76/fk5 eop constants
%    ttt         - julian centuries of tt                        centuries
%
%  outputs       :
%    rmod       - position vector mod                            km
%    vmod       - velocity vector mod                            km/s
%
%  locals        :
%    prec        - matrix for mod - eci
%
%  coupling      :
%   precess      - rotation for precession
%
%  references    :
%    vallado       2022, 226
%
% [reci,veci,aeci] = mod2eci  ( rmod,vmod,amod,ttt );
% ----------------------------------------------------------------------------

function [reci,veci,aeci] = mod2eci  ( rmod,vmod,amod,ttt )

    [prec,psia,wa,ea,xa] = precess ( ttt, '80' );

    reci = prec*rmod;

    veci = prec*vmod;

    aeci = prec*amod;

end