% ------------------------------------------------------------------------------
%
%                           function feci2NTW
%
%  defines coordinate rotation for ECI system to NTW
%  author        : sal alfano      719-573-2600   26 aug 2011
%
%  revisions
%                -
%
%  inputs          description
%    ra0           - ECI position vector (m)
%    va0           - ECI velocity vector (m)
%
%  outputs       :
%    rot_dv      - 3x3 rotation matrix from NTW to ECI
%
%  locals        :
%
%  rot_dv = feci2NTW(ra0 , a0);
% ------------------------------------------------------------------------------

function rot_dv = feci2NTW(ra0 , va0)

    % -------------------------  implementation   -----------------
    yprime = va0/norm(va0);
    h = cross(ra0,va0);
    zprime = h/norm(h);
    xprime = cross(yprime,zprime);
    rot_dv = horzcat(xprime,yprime,zprime)';

end