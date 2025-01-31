% ----------------------------------------------------------------------------
%
%                           function eci_teme
%
%  this function transforms a vector from the true equator mean equinox system,
%  (teme) to the mean equator mean equinox (j2000) system.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    rteme       - position vector of date
%                  true equator, mean equinox                     km
%    vteme       - velocity vector of date
%                  true equator, mean equinox                     km / s
%    iau80arr     - iau80 array of values
%    ttt         - julian centuries of tt                         centuries
%    ddpsi       - delta psi correction to gcrf                   rad
%    ddeps       - delta eps correction to gcrf                   rad
%
%  outputs       :
%    reci        - position vector eci                            km
%    veci        - velocity vector eci                            km / s
%
%  locals :
%    prec        - matrix for eci - mod
%    nutteme     - matrix for mod - teme - an approximation for nutation
%    eqeg        - rotation for equation of equinoxes(geometric terms only)
%    tm          - combined matrix for teme2eci
%
%  coupling      :
%   precess      - rotation for precession                        eci - mod
%   nutation     - rotation for nutation                          eci - tod
%
%  references :
%    vallado       2022, 233
%
% [rteme, vteme, ateme] = eci2teme  ( reci, veci, aeci, iau80arr, ttt, ddpsi, ddeps);
% ----------------------------------------------------------------------------

function [rteme, vteme, ateme] = eci2teme  ( reci, veci, aeci, iau80arr, ttt, ddpsi, ddeps)

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

    tm = eqe * nut' * prec';

    rteme = tm * reci;
    vteme = tm * veci;
    ateme = tm * aeci;

end