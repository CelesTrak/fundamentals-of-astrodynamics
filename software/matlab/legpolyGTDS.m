% ----------------------------------------------------------------------------
%
%                           function legpolyGTDS
%
%   this function finds the Legendre polynomials for the gravity field. fortran
%   and matlab implementations will have indicies of +1 as they start at 1. note
%   these are valid for normalized and unnormalized expressions, as long as
%   the order is less than about 85 (L + m). this is the GTDS approach.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs        description                                   range / units
%    latgc       - Geocentric lat of satellite                   pi to pi rad
%    degree      - degree of gravity field                        1..85
%    order       - order of gravity field                         1..85
%
%  outputs       :
%    LegArrGU    - array of unnormalized Legendre polynomials gtds
%
%  locals :
%    L,m         - degree and order indices
%
%  coupling      :
%   none
%
%  references :
%    vallado       2022, 600
%
%   [legarrGU] = legpolyGTDS (latgc, normArr, degree, order);
%  ----------------------------------------------------------------------------

function [legarrGU] = legpolyGTDS (latgc, degree, order)

    legarrGU = zeros(degree+3, order+3);

    sinlat = sin(latgc);

    % ---------------------------  GTDS approach --------------------------
    legarrGU(0+1, 0+1) = 1.0;
    legarrGU(0+1, 1+1) = 0.0;
    legarrGU(1+1, 0+1) = sinlat;
    legarrGU(1+1, 1+1) = cos(latgc);
    for L = 2 : degree
        Li = L + 1;
        for m = 0 : L
            mi = m + 1;
            legarrGU(Li, mi) = 0.0;
        end
        for m = 0 : L
            mi = m + 1;
            % Legendre functions, zonal gtds
            if (m == 0)
                legarrGU(Li, mi) = ((2.0*L-1) * legarrGU(2, 1) * legarrGU(Li-1, mi) - (L-1) * legarrGU(Li-2, mi)) / L;
            else
                % associated Legendre functions
                if m == L
                    legarrGU(Li, mi) = (2.0*L-1) * legarrGU(2, 2) * legarrGU(Li-1, mi-1);  % sectoral part
                else
                    legarrGU(Li, mi) = legarrGU(Li-2, mi) + (2.0*L-1) * legarrGU(2, 2) * legarrGU(Li-1, mi-1);  % tesseral part
                end
            end

        end   % for m
    end   % for L

end
