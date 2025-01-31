% ------------------------------------------------------------------------------
%
%                           function cubicspl
%
%  this function performs cubic splining of an input zero crossing
%  function in order to find function values.
%
%  author        : david vallado                  719-573-2600     2 feb 2004
%
%  revisions
%                -
%  inputs          description                    range / units
%    p1,p2,p3,p4 - function values used for splining
%    t1,t2,t3,t4 - time values used for splining
%
%  outputs       :
%    acu0..acu3  - splined polynomial coefficients. acu3 t^3, etc
%
%  locals        : none
%
%  coupling      :
%
%  references    :
%    vallado       2013, 559, 1034
%
% [acu0, acu1, acu2, acu3] = cubicspl(p1, p2, p3, p4);
% ------------------------------------------------------------------------------

function [acu0, acu1, acu2, acu3] = cubicspl(p1, p2, p3, p4)

    acu0 = p2;
    acu1 = -p1/3.0 - 0.5*p2 + p3 -p4/6.0;
    acu2 = 0.5*p1 - p2 + 0.5*p3;
    acu3 = -p1/6.0 + 0.5*p2 - 0.5*p3 + p4/6.0;

end
