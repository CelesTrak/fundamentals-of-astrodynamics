% ----------------------------------------------------------------------------
%
%                           function teme2eci
%
%  this function transforms a vector from the true equator mean equinox system,
%    (teme) to the mean equator mean equinox (j2000) system.
%
%  author        : david vallado                  719-573-2600   30 oct 2017
%
%  inputs          description                    range / units
%    rteme       - position vector of date
%                    true equator, mean equinox   km
%    vteme       - velocity vector of date
%                    true equator, mean equinox   km/s
%    ateme       - acceleration vector of date (if available, else set to 0)
%                    true equator, mean equinox   km/s2
%    ttt         - julian centuries of tt         centuries
%    ddpsi       - delta psi correction to gcrf   rad
%    ddeps       - delta eps correction to gcrf   rad
%
%  outputs       :
%    reci        - position vector eci            km
%    veci        - velocity vector eci            km/s
%    aeci        - acceleration vector eci        km/s2
%
%  locals        :
%    prec        - matrix for eci - mod
%    nutteme     - matrix for mod - teme - an approximation for nutation
%    eqeg        - rotation for equation of equinoxes (geometric terms only)
%    tm          - combined matrix for teme2eci
%
%  coupling      :
%   precess      - rotation for precession        eci - mod
%   nutation     - rotation for nutation          eci - tod
%
%  references    :
%    vallado       2013, 231-233
%
% [reci, veci, aeci] = teme2eci  ( rteme, vteme, ateme, iau80arr, ttt, ddpsi, ddeps);
% ----------------------------------------------------------------------------

function [reci, veci, aeci] = teme2eci  ( rteme, vteme, ateme, iau80arr, ttt, ddpsi, ddeps)

    [fArgs] = fundarg(ttt, '80');

    [prec,psia,wa,ea,xa] = precess ( ttt, '80' );

    [deltapsi, trueeps, meaneps, nut] = nutation  (ttt, ddpsi, ddeps, iau80arr, fArgs);

    % ------------------------ find eqeg ----------------------
    % rotate teme through just geometric terms
    eqeg = deltapsi* cos(meaneps);

    eqeg = rem (eqeg, 2.0*pi);

    eqe(1,1) =  cos(eqeg);
    eqe(1,2) =  sin(eqeg);
    eqe(1,3) =  0.0;
    eqe(2,1) = -sin(eqeg);
    eqe(2,2) =  cos(eqeg);
    eqe(2,3) =  0.0;
    eqe(3,1) =  0.0;
    eqe(3,2) =  0.0;
    eqe(3,3) =  1.0;

    tm = prec * nut * eqe';

    reci = tm * rteme;
    veci = tm * vteme;
    aeci = tm * ateme;

