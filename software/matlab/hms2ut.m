% ------------------------------------------------------------------------------
%
%                               procedure hms2ut
%
%  this procedure converts hours, minutes and seconds into universal time.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    hr          - hours                          0 .. 24
%    minute      - minutes                        0 .. 59
%    second      - seconds                        0.0 .. 59.99
%    direction   - which set of vars to output    from  too
%
%  outputs       :
%    ut          - universal time                 hrmin.sec
%
%  locals        :
%    none.
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2013, 199, alg 21, ex 3-10
%
% [ut] = hms2ut (hr, minute, second)
% -----------------------------------------------------------------------------

function [ut] = hms2ut (hr, minute, second)

    % ------------------------  implementation   ------------------
    ut = hr * 100.0 + minute + second* 0.01;
end  % hms2ut