% ----------------------------------------------------------------------------
%
%                           function trigpolyMont
%
%   this function finds the accumlated legendre polynomials and trigonmetric terms.
%   matlab/for implementations will have indicies of +1 because they start at 1.
%   beware that the multipliers must retain the original indices values!!
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs        description                                   range / units
%    recef       - satellite position vector, earth fixed        km
%    latgc       - latitude of satellite                         rad
%    degree      - size of gravity field                         1..2160..
%
%  outputs       :
%    VArr        - array of trig terms
%    WArr        - array of trig terms
%
%  locals :
%    L,m         - degree and order indices
%
%  coupling      :
%   none
%
%  references :
%    vallado       2022, 601
%
%  [VArr, WArr] = trigpolyMont(recef, latgc, degree);
% ----------------------------------------------------------------------------

function [VArr, WArr] = trigpolyMont(recef, latgc, degree)
    constastro;

    magr = mag(recef);

    % -------------------- montenbruck approach
    % now form first set of recursions for l=m on V and W
    % initial zonal values
    temp = re / (magr * magr);
    VArr(0+1, 0+1) = re / magr;
    %VArr(1+1, 0+1) = VArr(0+1, 0+1) * VArr(0+1, 0+1) * sin(latgc);
    VArr(1+1, 0+1) = VArr(0+1, 0+1) * VArr(0+1, 0+1) * recef(3) / magr;
    WArr(0+1, 0+1) = 0.0;
    WArr(1+1, 0+1) = 0.0;

    % find zonal terms
    for L = 2: degree+1
        Li = L + 1;
        mi = 0 + 1;
        VArr(Li, mi) = ((2*L-1) / L) * recef(3) * temp * VArr(Li-1, 0+1) - ((L - 1) / L) * temp * re * VArr(Li-1, 0+1);
        WArr(Li, mi) = 0.0;
    end

    % now tesseral and sectoral
    for L = 1: degree+1
        Li = L + 1;
        m = L;
        mi = m + 1;
        VArr(Li, mi) = (2*m-1) * recef(1) * temp * VArr(mi-1, mi-1) - recef(2) * temp * WArr(mi-1, mi-1);
        WArr(Li, mi) = (2*m-1) * recef(1) * temp * WArr(mi-1, mi-1) - recef(2) * temp * VArr(mi-1, mi-1);

        for (k = m + 1: degree + 1)
            ki = k + 1;
            VArr(ki, mi) = ((2*k-1) / (k-m)) * recef(3) * temp * VArr(ki-1, mi) - ...
                ((k+m-1) / (k-m)) * temp * re * VArr(ki-2, mi);
            WArr(ki, mi) = ((2*k-1) / (k-m)) * recef(3) * temp * WArr(ki-1, mi) - ...
                ((k+m-1) / (k-m)) * temp * re * WArr(ki-2, mi);
        end
    end

end
