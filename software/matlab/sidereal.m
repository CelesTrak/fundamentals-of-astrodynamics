% -----------------------------------------------------------------------------
%
%                           function sidereal
%
%  this function calculates the transformation matrix that accounts for the
%    effects of sidereal time. Notice that deltaspi should not be moded to a
%    positive number because it is multiplied rather than used in a
%    trigonometric argument.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    jdut1       - julian centuries of ut1                           days
%    deltapsi    - nutation angle                                    rad
%    meaneps     - mean obliquity of the ecliptic                    rad
%    fArgs06       - fundamental arguments in an array
%    lod         - length of day                                     sec
%    eqeterms    - terms for ast calculation                         0,2
%    opt         - method option                               e00cio, e96, e80
%
%  outputs       :
%    sidereal    - transformation matrix for pef - tod
%
%  locals        :
%    gmst         - mean greenwich sidereal time                 0 to 2pi rad
%    ast         - apparent gmst                                 0 to 2pi rad
%
%  coupling      :
%
%  references    :
%    vallado       2022, 214, 225
%
% [st, stdot]  = sidereal (jdut1, deltapsi, meaneps, omega, lod, eqeterms, opt );
% ----------------------------------------------------------------------------

function [st,stdot]  = sidereal (jdut1, deltapsi, meaneps, omega, lod, eqeterms, opt )
    constastro;

    if (opt == '80')
        % ------------------------ find gmst --------------------------
        gmst= gstime( jdut1 );

        % ------------------------ find mean ast ----------------------
        % after 1997, kinematic terms apply as well as gemoetric in eqe
        if (jdut1 > 2450449.5 ) && (eqeterms > 0)
            ast= gmst + deltapsi* cos(meaneps) ...
                + 0.00264*pi /(3600*180)*sin(omega) ...
                + 0.000063*pi /(3600*180)*sin(2.0 *omega);
        else
            ast= gmst + deltapsi* cos(meaneps);
        end

        ast = rem (ast, 2.0*pi);

        %fprintf(1,'st gmst %11.8f ast %11.8f ome  %11.8f \n', gmst*180/pi, ast*180/pi, omegaearth*180/pi );
    else
        % julian centuries of ut1
        tut1d= jdut1 - 2451545.0;

        era = twopi * ( 0.7790572732640 + 1.00273781191135448 * tut1d );
        era = rem (era,twopi);

        %if iauhelp == 'y'
        fprintf(1,'era%11.7f  \n',era*180/pi );
        %  end;
        ast = era;
    end

    st(1,1) =  cos(ast);
    st(1,2) = -sin(ast);
    st(1,3) =  0.0;
    st(2,1) =  sin(ast);
    st(2,2) =  cos(ast);
    st(2,3) =  0.0;
    st(3,1) =  0.0;
    st(3,2) =  0.0;
    st(3,3) =  1.0;

    % compute sidereal time rate matrix
    thetasa    = earthrot * (1.0  - lod/86400.0 );
    omegaearth = thetasa;

    stdot(1,1) = -omegaearth * sin(ast);
    stdot(1,2) = -omegaearth * cos(ast);
    stdot(1,3) =  0.0;
    stdot(2,1) =  omegaearth * cos(ast);
    stdot(2,2) = -omegaearth * sin(ast);
    stdot(2,3) =  0.0;
    stdot(3,1) =  0.0;
    stdot(3,2) =  0.0;
    stdot(3,3) =  0.0;

end