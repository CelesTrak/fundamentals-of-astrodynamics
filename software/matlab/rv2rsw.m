% -----------------------------------------------------------------------------
%
%                           function rv2rsw
%
%  this function converts position and velocity vectors into radial, along-
%    track, and MathTimeLibr.cross-track coordinates. note that sometimes the middle vector
%    is called in-track.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r           - position vector                          km
%    v           - velocity vector                          km/s
%
%  outputs       :
%    rrsw        - position vector                          km
%    vrsw        - velocity vector                          km/s
%    transmat    - transformation matrix
%
%  locals        :
%    tempvec     - temporary vector
%    rvec,svec,wvec - direction Math.Cosines
%
%  coupling      :
%
%  references    :
%    vallado       2022, 166
%
% [rrsw, vrsw,transmat] = rv2rsw( reci, veci );
% ------------------------------------------------------------------------------

function [rrsw, vrsw, transmat] = rv2rsw(reci, veci)

    % each of the components must be unit vectors
    % radial component
    rvec = unit(reci);

    % cross-track component
    wvec    = cross(reci, veci);
    wvec    = unit( wvec );

    % along-track component
    svec    = cross(wvec, rvec);
    svec    = unit( svec );

    % assemble transformation matrix from eci to rsw frame (individual
    % components arranged in row vectors)
    transmat(1,1) = rvec(1);
    transmat(1,2) = rvec(2);
    transmat(1,3) = rvec(3);
    transmat(2,1) = svec(1);
    transmat(2,2) = svec(2);
    transmat(2,3) = svec(3);
    transmat(3,1) = wvec(1);
    transmat(3,2) = wvec(2);
    transmat(3,3) = wvec(3);

    rrsw = matvecmult(transmat, reci, 3);
    vrsw = matvecmult(transmat, veci, 3);

    %   alt approach
    %       rrsw(1) = mag(reci);
    %       rrsw(2) = 0.0;
    %       rrsw(3) = 0.0;
    %       vrsw(1) = dot(reci,veci)/rrsw(1);
    %       vrsw(2) = sqrt(veci(1)^2 + veci(2)^2 + veci(3)^2 - vrsw(1)^2);
    %       vrsw(3) = 0.0;

end