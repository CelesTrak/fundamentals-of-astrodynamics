% ----------------------------------------------------------------------------
%
%                           function legpolyn
%
%   this function finds the Legendre polynomials for the gravity field. note that
%   the arrays are indexed from 0 to coincide with the usual nomenclature (eq 8-21
%   in my text). fortran and matlab implementations will have indicies of +1
%   as they start at 1. note that some recursions at high degree tesseral terms
%   experience error for resonant orbits. these are valid for normalized and
%   unnormalized expressions, as long as the remaining equations are consistent.
%   for satellite operations, orders up to about 120 are valid.
%
%  author        : david vallado                    719-573-2600  10 oct 2019
%
%  inputs        description                                   range / units
%    latgc       - Geocentric lat of satellite                   -pi to pi rad
%    lon         - longitude of satellite                        0 - 2pi rad
%    order       - size of gravity field                         1..~170
%
%  outputs       :
%    LegArrN     - array of normalized Legendre polynomials montenbruck
%    LegArrN1    - array of normalized Legendre polynomials gtds
%    trigArr     - array of trigonometric terms
%
%  locals :
%    L,m         - degree and order indices
%
%  coupling      :
%   none
%
%  references :
%    vallado       2013, 597, Eq 8-57
%
%   [P, legarrMU, legarrGU, legarrMN, legarrGN] = legpolyn (latgc, order);
%  ----------------------------------------------------------------------------*/

function [P, legarrMU, legarrGU, legarrMN, legarrGN] = legpolyn (latgc, order)
    %function [legarrMN, legarrGN] = legpolyn ( latgc, order)

    legarrGU = zeros(order+3,order+3);
    legarrGN = zeros(order+3,order+3);
    legarrMU = zeros(order+3,order+3);
    legarrMN = zeros(order+3,order+3);

    % increment indicies because Matlab always starts at 1!!
    % the +1 is shown to highlight where the indicies are different

    % -------------------- perform recursions ---------------------- }
    % --------------------  montenbruck approach
    legarrMU(0+1, 0+1) = 1.0;
    legarrMU(0+1, 1+1) = 0.0;
    legarrMU(1+1, 0+1) = sin(latgc);
    legarrMU(1+1, 1+1) = cos(latgc);

    % Legendre functions, zonal
    for L = 2: order
        Li = L + 1;
        legarrMU(Li, Li) = (2.0*L-1) * legarrMU(2, 2) * legarrMU(Li-1, Li-1);
    end

    % associated Legendre functions
    for L = 2: order
        Li = L + 1;
        for m = 0: L-1
            mi = m + 1;
            if L == m + 1
                legarrMU(Li, mi) = (2.0*m+1) * legarrMU(2, 1) * legarrMU(mi, mi);
            else
                legarrMU(Li, mi) = 1.0 / (L - m)  ...
                    * ((2.0*L-1) * legarrMU(2, 1) * legarrMU(Li-1, mi) - (L+m-1) * legarrMU(Li-2, mi));
            end
        end
    end

    % normalize the legendre polynomials
    %if order < 150
    for L = 0 : order
        Li = L + 1;
        for m = 0 : L
            mi = m + 1;
            % find normalized or unnormalized depending on which is already in file
            % note that above n = 170, the factorial will return 0, thus affecting the results!!!!
            if (m == 0)
                conv(Li, mi) = sqrt((factorial(L - m) * (2.0 * L + 1)) / factorial(L + m));
            else
                conv(Li, mi) = sqrt((factorial(L - m) * 2.0 * (2 * L + 1)) / factorial(L + m));
            end
            legarrMN(Li, mi) = conv(Li, mi) * legarrMU(Li, mi);

        end   % for m
    end   % for L
    %end


    % -------------------- gtds approach
    legarrGU(0+1, 0+1) = 1.0;
    legarrGU(0+1, 1+1) = 0.0;
    legarrGU(1+1, 0+1) = sin(latgc);
    legarrGU(1+1, 1+1) = cos(latgc);
    for L = 2 : order
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


    %if order < 150
    for L = 0 : order
        Li = L + 1;
        for m = 0 : L
            mi = m + 1;
            % find normalized or unnormalized depending on which is already in file
            % note that above n = 170, the factorial will return 0, thus affecting the results!!!!
            if (m == 0)
                conv1(Li, mi) = sqrt((factorial(L-m) * (2.0*L+1)) / factorial(L + m));
            else
                conv1(Li, mi) = sqrt((factorial(L-m) * 2.0 * (2*L+1)) / factorial(L + m));
            end
            legarrGN(Li, mi) = conv1(Li, mi) * legarrGU(Li, mi);

        end   % for m
    end   % for L
    %end

    % -------------------- try new normalization
    % use MIT example - which seems to give unnormalized answers????

    sinx = sin(latgc);
    cosx = cos(latgc);

    % now mit/etc approach for larger degree and order
    order = 1500;
    Pxn = zeros(order+3, order+3);
    nmax = order;
    %U = sin(latgc);  % R_F(3)/RMAG;
    Pxn(0+1, 0+1) = 1.0;
    % MITgcm has sin here (C/S later), Adamo too??, AMM Brazil has it backwards
    Pxn(1+1, 0+1) = sqrt(3.0) * sin(latgc);
    Pxn(1+1, 1+1) = sqrt(3.0) * cos(latgc);
    for L = 2:order
        Li = L + 1;

        for m = 0 : L
            mi = m + 1;
            if Li == mi
                Pxn(mi, mi) = cosx * sqrt((2*m+1)/(2*m)) * Pxn(mi-1,mi-1);
            else
                alpha = sqrt( (2*L+1)*(2*L-1) / ((L-m)*(L+m)) );
                beta = sqrt( (2*L+1)*(L-m-1)*(L+m-1) / ((2*L-3)*(L+m)*(L-m)) );
                Pxn(Li, mi) = alpha * sinx * Pxn(Li-1, mi) - beta * Pxn(Li-2, mi);
            end
        end
    end

    % CK SHum
    %how many ts? (how long is 3rd dimension of coeff matrix)
    % his approach is co-latitude
    t = pi*0.5 - latgc;
    k=length(t);

    % initialize P matrix
    P=zeros(order+1,order+1,k);

    % calculate u=sin(t)
    u=sin(t);

    % initial values
    % P_0,0
    P(1,1,1:k)=ones(1,k);
    % P_1,1
    P(2,2,1:k)=sqrt(3).*u;
    % P_1,0
    P(2,1,1:k)=sqrt(3).*cos(t);

    % sectorals = diagonals, where m=l
    for in=3:(order+1)
        % actual n
        n=in-1;
        % here and below, have to use reshape because funky 3rd dimension
        % matrix element-wise multiplication
        P(in,in,1:k)=sqrt((2.*n+1)./(2.*n)).*cos(t).*reshape(P(in-1,in-1,1:k),1,[]);

    end

    % iterate through l
    for in=3:(order+1)
        % actual l
        l=in-1;
        % for each l, iterate through m for m<l
        for jn=1:(in-1)
            % actual m
            m=jn-1;
            % dont have to save alpha and betas,
            % just use same variable for each iteration
            al=sqrt(((2.*l-1).*(2.*l+1))./((l-m).*(l+m)));
            be=sqrt(((2.*l+1).*(l+m-1).*(l-m-1))./((2.*l-3).*(l-m).*(l+m)));
            % compute P_l,m (remeber offset index)
            % again, funky reshapes
            P(in,jn,1:k)=al.*cos(t).*reshape(P(in-1,jn,1:k),1,[])-be.*reshape(P(in-2,jn,1:k),1,[]);
        end
    end

    % his approach is co-latitude
    t = pi*0.5 - latgc;
    n = order;
    m = order;
    %[p] = PBrockmeier(n,m,t);

end

function [p] = PBrockmeier(n,m,t)
    % Matlab Code [Credit: John Brockmeier, 2015]
    if (n == m) % P(n,n,t) :special cases
        if (n == 0) % P(0,0,t):n=0
            p = 1;
        elseif (n == 1) % P(1,1,t):n=1
            p = sqrt(3)*sqrt((1-(t^2)));
        else % P(n,n,t):n>=2
            p = sqrt(((2*n)+1)/(2*n))*sqrt((1-t^2))*PBrockmeier(n-1,n-1,t);
        end
    elseif (m == n-1) % P(n,n-1,t):n>=1
        p = sqrt((2*n)+1)*t*PBrockmeier(n-1,n-1,t);
    else % P(n,m,t) :0<=m<=n-2, n>=2
        p = A(n,m)*t*PBrockmeier(n-1,m,t)-B(n,m)*PBrockmeier(n-2,m,t);
    end
end

% Recursive 'constants' alpha & beta
function [a] = A(n,m)
    a = sqrt((((2*n)-1)*((2*n)+1))/((n-m)*(n+m)));
end

function [b] = B(n,m)
    b = sqrt((((2*n)+1)*(n+m-1)*(n-m-1))/(((2*n)-3)*(n-m)*(n+m)));
end



