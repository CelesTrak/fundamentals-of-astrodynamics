% -----------------------------------------------------------------------------
%
%                           function nutationqmod
%
%  this function calculates the transformation matrix that accounts for the
%    effects of nutation within the qmod paradigm. There are several assumptions
%    mentioned in the comments below.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    ttt         - julian centuries of tt
%    iau80arr    - array of iau80 values
%    fArgs06       - fundamental arguments in an array
%
%  outputs       :
%    deltapsi    - nutation in longiotude angle                   rad
%    trueeps     - true obliquity of the ecliptic                 rad
%    meaneps     - mean obliquity of the ecliptic                 rad
%    nut         - transform matrix for tod - mod
%
%  locals        :
%    iar80       - integers for fk5 1980
%    rar80       - reals for fk5 1980                             rad
%    fArgs06[0]    - delaunay element                               rad
%    fArgs06[1]    - delaunay element                               rad
%    fArgs06[2]    - delaunay element                               rad
%    fArgs06[3]    - delaunay element                               rad
%    fArgs06[4]    - delaunay element                               rad
%    deltaeps    - change in obliquity                            rad
%
%  coupling      :
%
%  references    :
%    vallado       2022, 224
%
% [deltapsi, trueeps, meaneps, omega,nut] = nutationqmod  (iau80arr, ttt )
% ----------------------------------------------------------------------------

function [deltapsi, trueeps, meaneps, omega, nut] = nutationqmod(iau80arr, ttt)
    deg2rad = pi/180.0;
    convrt= 0.0001 /3600.0;  % " to deg

    % ---- determine coefficients for iau 1980 nutation theory ----
    ttt2= ttt*ttt;
    ttt3= ttt2*ttt;

    % meaneps = -46.8150 *ttt - 0.00059 *ttt2 + 0.001813 *ttt3 + 84381.448;
    % reduced no extra terms
    meaneps = 84381.448;
    meaneps = rem( meaneps/3600.0 ,360.0  );
    meaneps = meaneps * deg2rad;

    [ fArgs96] = fundarg( ttt, '96' );
    omega = fArgs96(5);
    %fprintf(1,'nut del arg %11.7f  %11.7f  %11.7f  %11.7f  %11.7f  \n',l*180/pi,l1*180/pi,f*180/pi,d*180/pi,omega*180/pi );

    deltapsi= 0.0;
    deltaeps= 0.0;
    for i= 106:-1: 1  % not 106, but do they do it in reverse?? shouldn't be a big diff because only 9 terms
        tempval = iau80arr.iar80(i, 1) * fArgs96(1) + iau80arr.iar80(i, 2) * fArgs96(2) + iau80arr.iar80(i, 3) * fArgs96(3) + ...
            iau80arr.iar80(i, 4) * fArgs96(4) + iau80arr.iar80(i, 5) * fArgs96(5);
        deltapsi = deltapsi + (iau80arr.rar80(ii, 1) + iau80arr.rar80(ii, 2) * ttt) * sin(tempval);
        deltaeps = deltaeps + (iau80arr.rar80(ii, 3) + iau80arr.rar80(ii, 4) * ttt) * cos(tempval);
    end

    % --------------- find nutation parameters --------------------
    deltapsi = rem( deltapsi, 360.0  ) * deg2rad;  % + ddpsi/deg2rad,360.0
    deltaeps = rem( deltaeps, 360.0  ) * deg2rad;  % + ddeps/deg2rad,360.0
    trueeps  = meaneps + deltaeps;

    %fprintf(1,'meaneps %11.7f dp  %11.7f de  %11.7f te  %11.7f  \n',meaneps*180/pi,deltapsi*180/pi,deltaeps*180/pi,trueeps*180/pi );

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
    nut;

    %deltaeps =   3.543997691923817e-005

    %deltapsi =   -5.978331920759111e-005

    % correct way
    %   n1 = rot1mat( trueeps );
    %   n2 = rot3mat( deltapsi );
    %   n3 = rot1mat( -meaneps );
    %   nut = n3*n2*n1

    % eutelsat approximation
    n1 = rot1mat( deltaeps );
    n2 = rot2mat( -deltapsi * sineps );
    %        n3 = rot1mat( -meaneps );
    % or nut n2 n1??????????????
    nut = n2*n1;  % no n3...

    % simulate c# where the matrices are inverted already, but then mutlipled
    % in n2 n1 order.
    %n1 = rot1mat( -deltaeps );
    %n2 = rot2mat( deltapsi * sineps );
    %nut = n2'*n1'  % no n3...

    %pause;
    %%deltaeps
    %deltapsi

    %nut =[...
    %   0.999999998212977  -0.000054849510209  -0.000023781851926; ...
    %   0.000054850353003   0.999999997867747   0.000035439324690; ...
    %   0.000023779908046  -0.000035440629070   0.999999999089239];

    %         nut

end