% -----------------------------------------------------------------------------
%
%                           function iau06gst
%
%  this function finds the greenwich sidereal time (iau-2006/2000).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    jdut1       - julian date in ut1                       days from 4713 bc
%    deltapsi    - nutation angle                           rad
%    ttt         - julian centuries
%    iau06arr    - array of iau06 values
%    fArgs       - fundamental arguments                    rad
%
%  outputs       :
%    gst         - greenwich sidereal time                  0 to 2pi rad
%
%  coupling      :
%    none
%
%  references    :
%    vallado       2022, 217, eq 3-71
%
% [gst, st] = iau06gst(jdut1, deltapsi, ttt, fArgs06);
% -----------------------------------------------------------------------------

function [gst, st] = iau06gst(jdut1, deltapsi, ttt, fArgs06)
    sethelp;
    constastro;

    deg2rad = pi/180.0;
    % " to rad
    convrt  = pi / (180.0*3600.0);

    ttt2 = ttt  * ttt;
    ttt3 = ttt2 * ttt;
    ttt4 = ttt2 * ttt2;
    ttt5 = ttt3 * ttt2;

    % mean obliquity of the ecliptic
    epsa = 84381.406- 46.836769 * ttt - 0.0001831 * ttt2 + 0.00200340 * ttt3 - 0.000000576 * ttt4 ...
        - 0.0000000434 * ttt5; % "
    epsa = rem(epsa/3600.0 ,360.0  ); % deg
    epsa = epsa * deg2rad; % rad

    %  evaluate the ee complementary terms
    gstsum0 = 0.0;
    for i = 33: -1 : 1
        tempval = iau06arr.ag0i(i, 1) * fArgs06(1) + iau06arr.ag0i(i, 2) * fArgs06(2) + iau06arr.ag0i(i, 3) * fArgs06(3) ...
        + iau06arr.ag0i(i, 4) * fArgs06(4) + iau06arr.ag0i(i, 5) * fArgs06(5) + iau06arr.ag0i(i, 6) * fArgs06(6) ...
        + iau06arr.ag0i(i, 7) * fArgs06(7) + iau06arr.ag0i(i, 8) * fArgs06(8) + iau06arr.ag0i(i, 9) * fArgs06(9) ...
        + iau06arr.ag0i(i, 10) * fArgs06(10) + iau06arr.ag0i(i, 11) * fArgs06(11) + iau06arr.ag0i(i, 12) * fArgs06(12) ...
        + iau06arr.ag0i(i, 13) * fArgs06(13) + iau06arr.ag0i(i, 14) * fArgs06(14);
        gstsum0 = gstsum0 + iau06arr.ago(i,1)*sin(tempval) + iau06arr.ago(i,2)*cos(tempval); % rad
    end

    gstsum1 = 0.0;
    for j = 1: -1 : 1
        i = 33 + j;
        tempval = iau06arr.ag0i(i, 1) * fArgs06(1) + iau06arr.ag0i(i, 2) * fArgs06(2) + iau06arr.ag0i(i, 3) * fArgs06(3) ...
        + iau06arr.ag0i(i, 4) * fArgs06(4) + iau06arr.ag0i(i, 5) * fArgs06(5) + iau06arr.ag0i(i, 6) * fArgs06(6) ...
        + iau06arr.ag0i(i, 7) * fArgs06(7) + iau06arr.ag0i(i, 8) * fArgs06(8) + iau06arr.ag0i(i, 9) * fArgs06(9) ...
        + iau06arr.ag0i(i, 10) * fArgs06(10) + iau06arr.ag0i(i, 11) * fArgs06(11) + iau06arr.ag0i(i, 12) * fArgs06(12) ...
        + iau06arr.ag0i(i, 13) * fArgs06(13) + iau06arr.ag0i(i, 14) * fArgs06(14);
        gstsum1 = gstsum1 + iau06arr.ago(i,1)*ttt*sin(tempval) + iau06arr.ago(i,2)*ttt*cos(tempval);
    end

    eect2000 = gstsum0 + gstsum1 * ttt;  % rad

    % equation of the equinoxes
    ee2000 = deltapsi * cos(epsa) + eect2000;  % rad

    %  earth rotation angle
    tut1d= jdut1 - 2451545.0;
    era = twopi * ( 0.7790572732640 + 1.00273781191135448 * tut1d );
    era = rem (era,twopi);  % rad

    %  greenwich mean sidereal time, iau 2000.
    gmst2000 = era + (0.014506 + 4612.156534 * ttt + 1.3915817 * ttt2 ...
        - 0.00000044 * ttt3 + 0.000029956 * ttt4 + 0.0000000368 * ttt5) * convrt; % " to rad

    gst = gmst2000 + ee2000; % rad

    if iauhelp == 'y'
        fprintf(1,'meanobl %11.7f getsum %11.7f %11.7f eect %11.7f  \n',epsa*180/pi,gstsum0*180/pi,gstsum1*180/pi,eect2000*180/pi );
        fprintf(1,'ee2000 %11.7f gmst2000 %11.7f gst %11.7f  \n',ee2000*180/pi,gmst2000*180/pi,gst*180/pi );
    end

    % transformation matrix
    st(1,1) =  cos(gst);
    st(1,2) = -sin(gst);
    st(1,3) =  0.0;

    st(2,1) =  sin(gst);
    st(2,2) =  cos(gst);
    st(2,3) =  0.0;

    st(3,1) =  0.0;
    st(3,2) =  0.0;
    st(3,3) =  1.0;

end



