% -----------------------------------------------------------------------------
%
%                           function nutation
%
%  this function calculates the transformation matrix that accounts for the
%    effects of nutation within iau-76/fk5.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    ttt         - julian centuries of tt
%    ddpsi       - delta psi correction to gcrf                      rad
%    ddeps       - delta eps correction to gcrf                      rad
%    iau80arr    - array of iau80 values
%    fArgs06       - fundamental arguments in an array
%
%  outputs       :
%    deltapsi    - nutation in longitude angle                       rad
%    deltaeps    - nutation in obliquity angle                       rad
%    trueeps     - true obliquity of the ecliptic                    rad
%    meaneps     - mean obliquity of the ecliptic                    rad
%    nutation    - transform matrix for tod - mod
%
%  coupling      :
%
%  references    :
%    vallado       2022, 214, 225
%
% [deltapsi, trueeps, meaneps, nut] = nutation (ttt, ddpsi, ddeps, iau80arr, fArgs);
% ----------------------------------------------------------------------------

function [deltapsi, trueeps, meaneps, nut] = nutation  (ttt, ddpsi, ddeps, iau80arr, fArgs)
    deg2rad = pi/180.0;

    % ---- determine coefficients for iau 1980 nutation theory ----
    ttt2= ttt*ttt;
    ttt3= ttt2*ttt;

    meaneps = -46.8150 *ttt - 0.00059 *ttt2 + 0.001813 *ttt3 + 84381.448;
    meaneps = rem( meaneps/3600.0, 360.0 );
    meaneps = meaneps * deg2rad;

    deltapsi= 0.0;
    deltaeps= 0.0;
    for ii= 106:-1: 1
        tempval = iau80arr.iar80(ii, 1) * fArgs(1) + iau80arr.iar80(ii, 2) * fArgs(2) + iau80arr.iar80(ii, 3) * fArgs(3) + ...
            iau80arr.iar80(ii, 4) * fArgs(4) + iau80arr.iar80(ii, 5) * fArgs(5);
        deltapsi = deltapsi + (iau80arr.rar80(ii, 1) + iau80arr.rar80(ii, 2) * ttt) * sin(tempval);
        deltaeps = deltaeps + (iau80arr.rar80(ii, 3) + iau80arr.rar80(ii, 4) * ttt) * cos(tempval);
    end

    % --------------- find nutation parameters --------------------
    deltapsi = rem( deltapsi + ddpsi, 2.0 * pi );
    deltaeps = rem( deltaeps + ddeps, 2.0 * pi );
    trueeps  = meaneps + deltaeps;

    %fprintf(1,'meaneps %11.7f dp  %11.7f de  %11.7f te  %11.7f  ttt  %11.7f \n',meaneps*180/pi,deltapsi*180/pi,deltaeps*180/pi,trueeps*180/pi, ttt );

    cospsi  = cos(deltapsi);
    sinpsi  = sin(deltapsi);
    coseps  = cos(meaneps);
    sineps  = sin(meaneps);
    costrueeps = cos(trueeps);
    sintrueeps = sin(trueeps);

    nut(1,1) =  cospsi;
    nut(1,2) =  costrueeps * sinpsi;
    nut(1,3) =  sintrueeps * sinpsi;
    nut(2,1) = -coseps * sinpsi;
    nut(2,2) =  costrueeps * coseps * cospsi + sintrueeps * sineps;
    nut(2,3) =  sintrueeps * coseps * cospsi - sineps * costrueeps;
    nut(3,1) = -sineps * sinpsi;
    nut(3,2) =  costrueeps * sineps * cospsi - sintrueeps * coseps;
    nut(3,3) =  sintrueeps * sineps * cospsi + costrueeps * coseps;

    %         fprintf(1,'nut matrix \n');
    %         nut
    %         fprintf(1,'nut rotations \n');
    %         n1 = rot1mat( trueeps );
    %         n2 = rot3mat( deltapsi );
    %         n3 = rot1mat( -meaneps );
    %         nut1 = n3*n2*n1

end