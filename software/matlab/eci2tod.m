% ----------------------------------------------------------------------------
%
%                           function eci_tod
%
%  this function transforms between the eci mean equator mean equinox (gcrf), and
%    the true of date frame (tod).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    reci        - position vector eci                           km
%    veci        - velocity vector eci                           km/s
%    direct      - direction                                     eto, efrom
%    iau80arr    - iau76/fk5 eop constants
%    jdtt        - julian date of tt                             days from 4713 bc
%    jdftt       - fractional julian centuries of tt             days
%    jdut1       - julian date of ut1                            days from 4713 bc
%    lod         - excess length of day                          sec
%    ddpsi       - delta psi correction to eci                   rad
%    ddeps       - delta eps correction to eci                   rad
%
%  outputs       :
%    rtod       - position vector tod                            km
%    vtod       - velocity vector tod                            km/s
%
%  locals        :
%    deltapsi    - nutation angle                                rad
%    trueeps     - true obliquity of the ecliptic                rad
%    meaneps     - mean obliquity of the ecliptic                rad
%    prec        - matrix for mod - eci
%    nut         - matrix for tod - mod
%
%  coupling      :
%   precess      - rotation for precession
%   nutation     - rotation for nutation
%
%  references    :
%    vallado       2022, 225
%
% [rtod, vtod, atod] = eci2tod(reci, veci, aeci, iau80arr, ttt, ddpsi, ddeps, ddx, ddy )
% ----------------------------------------------------------------------------

function [rtod, vtod, atod] = eci2tod(reci, veci, aeci, iau80arr, ttt, ddpsi, ddeps)
    showit = 'n';

    [fArgs] = fundarg(ttt, '80');

    [prec,psia,wa,ea,xa] = precess ( ttt, '80' );

    [deltapsi, trueeps, meaneps, nut] = nutation  (ttt, ddpsi, ddeps, iau80arr, fArgs);

    if showit == 'y'
        conv = pi / (180.0*3600.0);
        fprintf(1,'dpsi %11.7f trueeps %11.7f mean eps %11.7f deltaeps %11.7f \n', deltapsi/conv, trueeps/conv, meaneps/conv, (trueeps-meaneps)/conv);
        fprintf(1,'nut iau 76 \n');
        fprintf(1,'%20.14f %20.14f %20.14f \n',nut );
    end

    rtod=nut'*prec'*reci;

    vtod=nut'*prec'*veci;

    atod=nut'*prec'*aeci;

end