% -----------------------------------------------------------------------------
%
%                           function polarm
%
%  this function calculates the transformation matrix that accounts for polar
%    motion within the iau-76/fk5, iau-2000a, and iau2006/2000 equinox systems.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    xp          - polar motion coefficient                         rad
%    yp          - polar motion coefficient                         rad
%    ttt         - julian centuries of tt (00 theory only)
%    opt         - method option                                   e80, e96, e06cio
%
%  outputs       :
%    polarm      - transformation matrix for itrf - pef
%
%  locals        :
%    convrt      - conversion from arcsec to rad
%    sp          - s prime value (00 theory only)
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2022, 213, 224
%
% [pm] = polarm ( xp, yp, ttt, opt );
% ----------------------------------------------------------------------------

function [pm] = polarm ( xp, yp, ttt, opt )
    cosxp = cos(xp);
    sinxp = sin(xp);
    cosyp = cos(yp);
    sinyp = sin(yp);

    if (opt == '80')
        pm(1,1) =  cosxp;
        pm(1,2) =  0.0;
        pm(1,3) = -sinxp;
        pm(2,1) =  sinxp * sinyp;
        pm(2,2) =  cosyp;
        pm(2,3) =  cosxp * sinyp;
        pm(3,1) =  sinxp * cosyp;
        pm(3,2) = -sinyp;
        pm(3,3) =  cosxp * cosyp;

        % a1 = rot2mat(xp);
        % a2 = rot1mat(yp);
        % pm = a2*a1;
        % Approximate matrix using small angle approximations
        %pm(1,1) =  1.0;
        %pm(2,1) =  0.0;
        %pm(3,1) =  xp;
        %pm(1,2) =  0.0;
        %pm(2,2) =  1.0;
        %pm(3,2) = -yp;
        %pm(1,3) = -xp;
        %pm(2,3) =  yp;
        %pm(3,3) =  1.0;
    else
        convrt = pi / (3600.0*180.0);
        % approximate sp value in rad
        sp = -47.0e-6 * ttt * convrt;
        fprintf(1,'xp %16.14f, %16.14f sp %16.14g \n',xp, yp, sp);
        cossp = cos(sp);
        sinsp = sin(sp);

        %fprintf(1,' sp  %14.11f mas \n', sp/convrt);

        % form the matrix
        pm(1,1) =  cosxp * cossp;
        pm(1,2) = -cosyp * sinsp + sinyp * sinxp * cossp;
        pm(1,3) = -sinyp * sinsp - cosyp * sinxp * cossp;
        pm(2,1) =  cosxp * sinsp;
        pm(2,2) =  cosyp * cossp + sinyp * sinxp * sinsp;
        pm(2,3) =  sinyp * cossp - cosyp * sinxp * sinsp;
        pm(3,1) =  sinxp;
        pm(3,2) = -sinyp * cosxp;
        pm(3,3) =  cosyp * cosxp;

        % a1 = rot1mat(yp);
        % a2 = rot2mat(xp);
        % a3 = rot3mat(-sp);
        % pm = a3*a2*a1;
    end

end