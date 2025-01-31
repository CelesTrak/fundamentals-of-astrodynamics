% ------------------------------------------------------------------------------
%
%                               procedure ut2hms
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
% [hr, minute, second] = ut2hms (ut)
% -----------------------------------------------------------------------------

function [hr, minute, second] = ut2hms (ut)

    % ------------------------  implementation   ------------------
    hr = floor(ut * 0.01);
    minute = floor(ut - hr * 100.0);
    second= (ut - hr * 100.0 - minute) * 100.0;
end  % hms_ut