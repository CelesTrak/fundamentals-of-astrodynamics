% ----------------------------------------------------------------------------
%
%                           function getnormGottN
%
%   this function finds the normalization constant values for the Gottlieb approach.
%     these functions can be found once before any calls to the gottlieb method.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                                         range / units
%    degree      - size of gravity field (zonals)                          1.. 2000 ..
%
%  outputs       :
%    normxx   - arrays of normalized Legendre polynomials, constant portion
%
%  locals :
%    L,m         - degree and order indices
%
%  coupling      :
%   none
%
%  references :
%    eckman, brown, adamo 2016 paper and code in matlab
%    vallado       2022, 600, Eq 8-56
%
% [norm1, norm2, norm11, normn10, norm1m, norm2m, normn1] = getnormGottN(degree);
% ---------------------------------------------------------------------------

function [norm1, norm2, norm11, normn10, norm1m, norm2m, normn1] = getnormGottN(degree)

    % this can be done ahead of time
    % these are 0 based arrays since they start at 2
    for L = 2:degree+1 %RAE
        norm1(L) = sqrt((2*L+1.0) / (2*L-1.0)); % eq 3-1 RAE
        norm2(L) = sqrt((2*L+1.0) / (2*L-3.0)); % eq 3-2 RAE
        norm11(L) = sqrt((2*L+1.0) / (2*L))/(2*L-1.0); % eq 3-3 RAE
        normn10(L) = sqrt((L+1.0)*L * 0.5); % RAE

        for m = 1:L %RAE
            norm1m(L,m) = sqrt((L-m)*(2*L+1.0) / ((L+m)*(2*L-1.0))); % eq 3-4 RAE
            norm2m(L,m) = sqrt((L-m)*(L-m-1.0)*(2*L+1.0) / ((L+m)*(L+m-1.0)*(2*L-3.0))); % eq 3-5 RAE
            normn1(L,m) = sqrt((L+m+1.0)*(L-m)); % part of eq 3-9 RAE
        end %RAE
    end %RAE

end