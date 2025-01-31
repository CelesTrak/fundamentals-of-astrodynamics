% -----------------------------------------------------------------------------
%
%                           function cubicinterp
%
%  this function performs a cubic spline. four points are needed.
%
%  author        : david vallado                  719-573-2600   1 dec  2005
%
%  revisions
%
%  inputs          description                    range / units
%    valuein     - kp
%
%  outputs       :
%    out         - ap
%
%  locals        :
%                -
%
%  coupling      :
%    cubicspl
%
%  references    :
%    vallado       2013, 1027
%
% [answer] =  cubicinterp(p1a, p1b, p1c, p1d, p2a, p2b, p2c, p2d, valuein );
% --------------------------------------------------------------------------- 

function [answer] = cubicinterp(p1a, p1b, p1c, p1d, p2a, p2b, p2c, p2d, valuein )
    % -------- assign function points ---------
    [ac0, ac1, ac2, ac3] = cubicspl(p1a, p1b, p1c, p1d);
    [kc0, kc1, kc2, kc3] = cubicspl(p2a, p2b, p2c, p2d);

    % recover the original function values
    % use the normalized time first, but at an arbitrary interval
    [r1r, r1i, r2r, r2i, r3r, r3i] = cubic(kc3, kc2, kc1, kc0 - valuein, 'R');

    if ((r1r >= -0.000001) && (r1r <= 1.001))
        value = r1r;
    else
        if ((r2r >= -0.000001) && (r2r <= 1.001))
            value = r2r;
        else
            if ((r3r >= -0.000001) && (r3r <= 1.001))
                value = r3r;
            else
                value = 0.0;
                fprintf(1,'error in cubicinterp root %11.7f %11.7f %11.7f %11.7f \n', valuein, r1r, r2r, r3r);
            end
        end
    end

    answer = ac3 * value^3 + ac2 * value * value + ac1 * value + ac0;
    % cubicinterp

