% ---------------------------------------------------------------------------
%
%                           function kp2ap
%
%  this function converts ap to kp using cubic splines. notice the arrays go
%  beyond the range of values to permit endpoint evaluations without additional
%  logic. the arrays have an extra element so they will start at 1. also, the normal
%  cubic splining is usually between pts 2 and 3, but here it is between 3 and 4
%  because i've added 2 additional points at the start.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    apin        - ap
%
%  outputs       :
%    kpout       - kp
%
%  locals        :
%    idx         - index of function value above the input value so the input
%                  value is between the 2nd and 3rd point
%
%  coupling      :
%    cubicspl    - perform the splining operation given 4 points
%
%  references    :
%    vallado       2013, 558
%
% kpout = ap2kp(apin);
% ---------------------------------------------------------------------------

function kpout = ap2kp(apin)

    bap = [0 -0.00001 -0.001 0 2 3 4 5 6 7 9 12 15 18 22 27 32 39 48 56 67 80 94 111 132 154 179 207 236 300 400 900];
    bkp = [0 -6.6666667 -3.3333 0.0 3.3333 6.6667 10 13.3333 16.6667 20 23.3333 26.6667 ...
        30 33.3333 36.6667 40 43.3333 46.6667 50 53.3333 56.6667 ...
        60 63.3333 66.6667 70 73.3333 76.6667 80 83.3333 86.6667  90 93.3333];

    kpout = 0.0;

    % find starting point in files based on input
    idx = 1;
    while ((idx < 33) && (apin > bap(idx)))
        idx = idx + 1;
    end

    if idx > 2
        kpout = cubicinterp(bkp(idx - 2), bkp(idx - 1), bkp(idx), bkp(idx + 1), ...
            bap(idx - 2), bap(idx - 1), bap(idx), bap(idx + 1), ...
            apin );
    else
        kpout = 0.0;
    end

end
