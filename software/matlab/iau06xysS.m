% ----------------------------------------------------------------------------
%
%                           function iau06xysS
%
%  this function calculates the XYS parameters for the iau2006 cio theory directly
%  from the summations.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    ttt         - julian centuries of tt
%    iau06arr    - iau06 array of values
%    fArgs06       - fundamental arguments in an array
%    1  l        - mean anomaly of the moon                              rad
%    2  l1       - mean anomaly of the Sun                               rad
%    3  f        - mean longitude of the Moon minus that of asc node     rad
%    4  d        - mean elongation of the Moon from the Sun              rad
%    5  omega    - mean longitude of the ascending node of the Moon      rad
%    6-14  planetary longitudes                                          rad
%
%  outputs       :
%    x           - coordinate of cio                                     rad
%    y           - coordinate of cio                                     rad
%    s           - coordinate                                            rad
%
%  references    :
%    vallado       2022, 214-216
%
% [x ,y, s] = iau06xysS (iau06arr, fArgs06, ttt );
% ----------------------------------------------------------------------------

function [x ,y, s] = iau06xysS (iau06arr, fArgs06, ttt )
    sethelp;

    % " to rad
    convrt  = pi / (180.0*3600.0);
    deg2rad = pi / 180.0;

    ttt2 = ttt  * ttt;
    ttt3 = ttt2 * ttt;
    ttt4 = ttt2 * ttt2;
    ttt5 = ttt3 * ttt2;

    % ---------------- first find x
    % the iers code puts the constants in here, however
    % don't sum constants in here because they're larger than the last few terms
    xsum0 = 0.0;
    for i = 1306: -1 : 1
        tempval = fArgs06(1)*iau06arr.ax0i(i,1) + iau06arr.ax0i(i,2)*fArgs06(2) + iau06arr.ax0i(i,3)*fArgs06(3) + iau06arr.ax0i(i,4)*fArgs06(4) + iau06arr.ax0i(i,5)*fArgs06(5) + ...
            iau06arr.ax0i(i,6)*fArgs06(6)  + iau06arr.ax0i(i,7)*fArgs06(7)  + iau06arr.ax0i(i,8)*fArgs06(8)  + iau06arr.ax0i(i,9)*fArgs06(9) + ...
            iau06arr.ax0i(i,10)*fArgs06(10) + iau06arr.ax0i(i,11)*fArgs06(11) + iau06arr.ax0i(i,12)*fArgs06(12) + iau06arr.ax0i(i,13)*fArgs06(13) + iau06arr.ax0i(i,14)*fArgs06(14);
        xsum0 = xsum0 + iau06arr.ax0(i,1) * sin(tempval) + iau06arr.ax0(i,2) * cos(tempval);
    end
    xsum1 = 0.0;
    % note that the index changes here to j. this is because the iau06arr.ax0i etc
    % indicies go from 1 to 1600, but there are 5 groups. the i index counts through each
    % calculation, and j takes care of the individual summations. note that
    % this same process is used for y and s.
    for j = 253: -1 : 1
        i = 1306 + j;
        tempval = fArgs06(1)*iau06arr.ax0i(i,1) + iau06arr.ax0i(i,2)*fArgs06(2) + iau06arr.ax0i(i,3)*fArgs06(3) + iau06arr.ax0i(i,4)*fArgs06(4) + iau06arr.ax0i(i,5)*fArgs06(5) + ...
            iau06arr.ax0i(i,6)*fArgs06(6)  + iau06arr.ax0i(i,7)*fArgs06(7) + iau06arr.ax0i(i,8)*fArgs06(8)  + iau06arr.ax0i(i,9)*fArgs06(9) + ...
            iau06arr.ax0i(i,10)*fArgs06(10) + iau06arr.ax0i(i,11)*fArgs06(11) + iau06arr.ax0i(i,12)*fArgs06(12) + iau06arr.ax0i(i,13)*fArgs06(13) + iau06arr.ax0i(i,14)*fArgs06(14);
        xsum1 = xsum1 + iau06arr.ax0(i,1)*sin(tempval) + iau06arr.ax0(i,2)*cos(tempval);
    end
    xsum2 = 0.0;
    for j = 36: -1 : 1
        i = 1306 + 253 + j;
        tempval = fArgs06(1)*iau06arr.ax0i(i,1) + iau06arr.ax0i(i,2)*fArgs06(2) + iau06arr.ax0i(i,3)*fArgs06(3) + iau06arr.ax0i(i,4)*fArgs06(4) + iau06arr.ax0i(i,5)*fArgs06(5) + ...
            iau06arr.ax0i(i,6)*fArgs06(6)  + iau06arr.ax0i(i,7)*fArgs06(7)  + iau06arr.ax0i(i,8)*fArgs06(8)  + iau06arr.ax0i(i,9)*fArgs06(9) + ...
            iau06arr.ax0i(i,10)*fArgs06(10) + iau06arr.ax0i(i,11)*fArgs06(11) + iau06arr.ax0i(i,12)*fArgs06(12) + iau06arr.ax0i(i,13)*fArgs06(13) + iau06arr.ax0i(i,14)*fArgs06(14);
        xsum2 = xsum2 + iau06arr.ax0(i,1)*sin(tempval) + iau06arr.ax0(i,2)*cos(tempval);
    end
    xsum3 = 0.0;
    for j = 4: -1 : 1
        i = 1306 + 253 + 36 + j;
        tempval = fArgs06(1)*iau06arr.ax0i(i,1) + iau06arr.ax0i(i,2)*fArgs06(2) + iau06arr.ax0i(i,3)*fArgs06(3) + iau06arr.ax0i(i,4)*fArgs06(4) + iau06arr.ax0i(i,5)*fArgs06(5) + ...
            iau06arr.ax0i(i,6)*fArgs06(6)  + iau06arr.ax0i(i,7)*fArgs06(7)  + iau06arr.ax0i(i,8)*fArgs06(8)  + iau06arr.ax0i(i,9)*fArgs06(9) + ...
            iau06arr.ax0i(i,10)*fArgs06(10) + iau06arr.ax0i(i,11)*fArgs06(11) + iau06arr.ax0i(i,12)*fArgs06(12) + iau06arr.ax0i(i,13)*fArgs06(13) + iau06arr.ax0i(i,14)*fArgs06(14);
        xsum3 = xsum3 + iau06arr.ax0(i,1)*sin(tempval) + iau06arr.ax0(i,2)*cos(tempval);
    end
    xsum4 = 0.0;
    for j = 1: -1 : 1
        i = 1306 + 253 + 36 + 4 + j;
        tempval = fArgs06(1)*iau06arr.ax0i(i,1) + iau06arr.ax0i(i,2)*fArgs06(2) + iau06arr.ax0i(i,3)*fArgs06(3) + iau06arr.ax0i(i,4)*fArgs06(4) + iau06arr.ax0i(i,5)*fArgs06(5) + ...
            iau06arr.ax0i(i,6)*fArgs06(6)  + iau06arr.ax0i(i,7)*fArgs06(7)  + iau06arr.ax0i(i,8)*fArgs06(8)  + iau06arr.ax0i(i,9)*fArgs06(9) + ...
            iau06arr.ax0i(i,10)*fArgs06(10) + iau06arr.ax0i(i,11)*fArgs06(11) + iau06arr.ax0i(i,12)*fArgs06(12) + iau06arr.ax0i(i,13)*fArgs06(13) + iau06arr.ax0i(i,14)*fArgs06(14);
        xsum4 = xsum4 + iau06arr.ax0(i,1)*sin(tempval) + iau06arr.ax0(i,2)*cos(tempval);
    end

    x = -0.016617 + 2004.191898 * ttt - 0.4297829 * ttt2 ...
        - 0.19861834 * ttt3 - 0.000007578 * ttt4 + 0.0000059285 * ttt5; % "
    x = x*convrt + xsum0 + xsum1*ttt + xsum2*ttt2 + xsum3*ttt3 + xsum4*ttt4;  % rad

    if iauhelp == 'y'
        fprintf(1,'x %14.12f  %14.12f  %14.12f  %14.12f  %14.12f \n',xsum0,xsum1,xsum2,xsum3,xsum4 );
    end

    % ---------------- now find y
    ysum0 = 0.0;
    for i = 962: -1 : 1
        tempval = fArgs06(1)*iau06arr.ay0i(i,1) + iau06arr.ay0i(i,2)*fArgs06(2) + iau06arr.ay0i(i,3)*fArgs06(3) + iau06arr.ay0i(i,4)*fArgs06(4) + iau06arr.ay0i(i,5)*fArgs06(5) + ...
            iau06arr.ay0i(i,6)*fArgs06(6)  + iau06arr.ay0i(i,7)*fArgs06(7)  + iau06arr.ay0i(i,8)*fArgs06(8)  + iau06arr.ay0i(i,9)*fArgs06(9) + ...
            iau06arr.ay0i(i,10)*fArgs06(10) + iau06arr.ay0i(i,11)*fArgs06(11) + iau06arr.ay0i(i,12)*fArgs06(12) + iau06arr.ay0i(i,13)*fArgs06(13) + iau06arr.ay0i(i,14)*fArgs06(14);
        ysum0 = ysum0 + iau06arr.ay0(i,1)*sin(tempval) + iau06arr.ay0(i,2) * cos(tempval);
    end

    ysum1 = 0.0;
    for j = 277: -1 : 1
        i = 962 + j;
        tempval = fArgs06(1)*iau06arr.ay0i(i,1) + iau06arr.ay0i(i,2)*fArgs06(2) + iau06arr.ay0i(i,3)*fArgs06(3) + iau06arr.ay0i(i,4)*fArgs06(4) + iau06arr.ay0i(i,5)*fArgs06(5) + ...
            iau06arr.ay0i(i,6)*fArgs06(6)  + iau06arr.ay0i(i,7)*fArgs06(7)  + iau06arr.ay0i(i,8)*fArgs06(8)  + iau06arr.ay0i(i,9)*fArgs06(9) + ...
            iau06arr.ay0i(i,10)*fArgs06(10) + iau06arr.ay0i(i,11)*fArgs06(11) + iau06arr.ay0i(i,12)*fArgs06(12) + iau06arr.ay0i(i,13)*fArgs06(13) + iau06arr.ay0i(i,14)*fArgs06(14);
        ysum1 = ysum1 + iau06arr.ay0(i,1)*sin(tempval) + iau06arr.ay0(i,2)*cos(tempval);
    end
    ysum2 = 0.0;
    for j = 30: -1 : 1
        i = 962 + 277 + j;
        tempval = fArgs06(1)*iau06arr.ay0i(i,1) + iau06arr.ay0i(i,2)*fArgs06(2) + iau06arr.ay0i(i,3)*fArgs06(3) + iau06arr.ay0i(i,4)*fArgs06(4) + iau06arr.ay0i(i,5)*fArgs06(5) + ...
            iau06arr.ay0i(i,6)*fArgs06(6)  + iau06arr.ay0i(i,7)*fArgs06(7)  + iau06arr.ay0i(i,8)*fArgs06(8)  + iau06arr.ay0i(i,9)*fArgs06(9) + ...
            iau06arr.ay0i(i,10)*fArgs06(10) + iau06arr.ay0i(i,11)*fArgs06(11) + iau06arr.ay0i(i,12)*fArgs06(12) + iau06arr.ay0i(i,13)*fArgs06(13) + iau06arr.ay0i(i,14)*fArgs06(14);
        ysum2 = ysum2 + iau06arr.ay0(i,1)*sin(tempval) + iau06arr.ay0(i,2)*cos(tempval);
    end
    ysum3 = 0.0;
    for j = 5: -1 : 1
        i = 962 + 277 + 30 + j;
        tempval = fArgs06(1)*iau06arr.ay0i(i,1) + iau06arr.ay0i(i,2)*fArgs06(2) + iau06arr.ay0i(i,3)*fArgs06(3) + iau06arr.ay0i(i,4)*fArgs06(4) + iau06arr.ay0i(i,5)*fArgs06(5) + ...
            iau06arr.ay0i(i,6)*fArgs06(6)  + iau06arr.ay0i(i,7)*fArgs06(7)  + iau06arr.ay0i(i,8)*fArgs06(8)  + iau06arr.ay0i(i,9)*fArgs06(9) + ...
            iau06arr.ay0i(i,10)*fArgs06(10) + iau06arr.ay0i(i,11)*fArgs06(11) + iau06arr.ay0i(i,12)*fArgs06(12) + iau06arr.ay0i(i,13)*fArgs06(13) + iau06arr.ay0i(i,14)*fArgs06(14);
        ysum3 = ysum3 + iau06arr.ay0(i,1)*sin(tempval) + iau06arr.ay0(i,2)*cos(tempval);
    end
    ysum4 = 0.0;
    for j = 1: -1 : 1
        i = 962 + 277 + 30 + 5 + j;
        tempval = fArgs06(1)*iau06arr.ay0i(i,1) + iau06arr.ay0i(i,2)*fArgs06(2) + iau06arr.ay0i(i,3)*fArgs06(3) + iau06arr.ay0i(i,4)*fArgs06(4) + iau06arr.ay0i(i,5)*fArgs06(5) + ...
            iau06arr.ay0i(i,6)*fArgs06(6)  + iau06arr.ay0i(i,7)*fArgs06(7)  + iau06arr.ay0i(i,8)*fArgs06(8)  + iau06arr.ay0i(i,9)*fArgs06(9) + ...
            iau06arr.ay0i(i,10)*fArgs06(10) + iau06arr.ay0i(i,11)*fArgs06(11) + iau06arr.ay0i(i,12)*fArgs06(12) + iau06arr.ay0i(i,13)*fArgs06(13) + iau06arr.ay0i(i,14)*fArgs06(14);
        ysum4 = ysum4 + iau06arr.ay0(i,1)*sin(tempval) + iau06arr.ay0(i,2)*cos(tempval);
    end

    y = -0.006951 - 0.025896 * ttt - 22.4072747 * ttt2 ...
        + 0.00190059 * ttt3 + 0.001112526 * ttt4 + 0.0000001358 * ttt5;
    y = y*convrt + ysum0 + ysum1*ttt + ysum2*ttt2 + ysum3*ttt3 + ysum4*ttt4;  % rad

    if iauhelp == 'y'
        fprintf(1,'y %14.12f  %14.12f  %14.12f  %14.12f  %14.12f \n',ysum0/deg2rad,ysum1/deg2rad,ysum2/deg2rad,ysum3/deg2rad,ysum4/deg2rad );
    end

    % ---------------- now find s
    ssum0 = 0.0;
    for i = 33: -1 : 1
        tempval = fArgs06(1)*iau06arr.as0i(i,1) + iau06arr.as0i(i,2)*fArgs06(2) + iau06arr.as0i(i,3)*fArgs06(3) + iau06arr.as0i(i,4)*fArgs06(4) + iau06arr.as0i(i,5)*fArgs06(5) + ...
            iau06arr.as0i(i,6)*fArgs06(6)  + iau06arr.as0i(i,7)*fArgs06(7)  + iau06arr.as0i(i,8)*fArgs06(8)  + iau06arr.as0i(i,9)*fArgs06(9) + ...
            iau06arr.as0i(i,10)*fArgs06(10) + iau06arr.as0i(i,11)*fArgs06(11) + iau06arr.as0i(i,12)*fArgs06(12) + iau06arr.as0i(i,13)*fArgs06(13) + iau06arr.as0i(i,14)*fArgs06(14);
        ssum0 = ssum0 + iau06arr.as0(i,1)*sin(tempval) + iau06arr.as0(i,2)*cos(tempval);
    end
    ssum1 = 0.0;
    for j = 3: -1 : 1
        i = 33 + j;
        tempval = fArgs06(1)*iau06arr.as0i(i,1) + iau06arr.as0i(i,2)*fArgs06(2) + iau06arr.as0i(i,3)*fArgs06(3) + iau06arr.as0i(i,4)*fArgs06(4) + iau06arr.as0i(i,5)*fArgs06(5) + ...
            iau06arr.as0i(i,6)*fArgs06(6)  + iau06arr.as0i(i,7)*fArgs06(7)  + iau06arr.as0i(i,8)*fArgs06(8)  + iau06arr.as0i(i,9)*fArgs06(9) + ...
            iau06arr.as0i(i,10)*fArgs06(10) + iau06arr.as0i(i,11)*fArgs06(11) + iau06arr.as0i(i,12)*fArgs06(12) + iau06arr.as0i(i,13)*fArgs06(13) + iau06arr.as0i(i,14)*fArgs06(14);
        ssum1 = ssum1 + iau06arr.as0(i,1)*sin(tempval) + iau06arr.as0(i,2)*cos(tempval);
    end
    ssum2 = 0.0;
    for j = 25: -1 : 1
        i = 33 + 3 + j;
        tempval = fArgs06(1)*iau06arr.as0i(i,1) + iau06arr.as0i(i,2)*fArgs06(2) + iau06arr.as0i(i,3)*fArgs06(3) + iau06arr.as0i(i,4)*fArgs06(4) + iau06arr.as0i(i,5)*fArgs06(5) + ...
            iau06arr.as0i(i,6)*fArgs06(6)  + iau06arr.as0i(i,7)*fArgs06(7)  + iau06arr.as0i(i,8)*fArgs06(8)  + iau06arr.as0i(i,9)*fArgs06(9) + ...
            iau06arr.as0i(i,10)*fArgs06(10) + iau06arr.as0i(i,11)*fArgs06(11) + iau06arr.as0i(i,12)*fArgs06(12) + iau06arr.as0i(i,13)*fArgs06(13) + iau06arr.as0i(i,14)*fArgs06(14);
        ssum2 = ssum2 + iau06arr.as0(i,1)*sin(tempval) + iau06arr.as0(i,2)*cos(tempval);
    end
    ssum3 = 0.0;
    for j = 4: -1 : 1
        i = 33 + 3 + 25 + j;
        tempval = fArgs06(1)*iau06arr.as0i(i,1) + iau06arr.as0i(i,2)*fArgs06(2) + iau06arr.as0i(i,3)*fArgs06(3) + iau06arr.as0i(i,4)*fArgs06(4) + iau06arr.as0i(i,5)*fArgs06(5) + ...
            iau06arr.as0i(i,6)*fArgs06(6)  + iau06arr.as0i(i,7)*fArgs06(7)  + iau06arr.as0i(i,8)*fArgs06(8)  + iau06arr.as0i(i,9)*fArgs06(9) + ...
            iau06arr.as0i(i,10)*fArgs06(10) + iau06arr.as0i(i,11)*fArgs06(11) + iau06arr.as0i(i,12)*fArgs06(12) + iau06arr.as0i(i,13)*fArgs06(13) + iau06arr.as0i(i,14)*fArgs06(14);
        ssum3 = ssum3 + iau06arr.as0(i,1)*sin(tempval) + iau06arr.as0(i,2)*cos(tempval);
    end
    ssum4 = 0.0;
    for j = 1: -1 : 1
        i = 33 + 3 + 25 + 4 + j;
        tempval = fArgs06(1)*iau06arr.as0i(i,1) + iau06arr.as0i(i,2)*fArgs06(2) + iau06arr.as0i(i,3)*fArgs06(3) + iau06arr.as0i(i,4)*fArgs06(4) + iau06arr.as0i(i,5)*fArgs06(5) + ...
            iau06arr.as0i(i,6)*fArgs06(6)  + iau06arr.as0i(i,7)*fArgs06(7)  + iau06arr.as0i(i,8)*fArgs06(8)  + iau06arr.as0i(i,9)*fArgs06(9) + ...
            iau06arr.as0i(i,10)*fArgs06(10) + iau06arr.as0i(i,11)*fArgs06(11) + iau06arr.as0i(i,12)*fArgs06(12) + iau06arr.as0i(i,13)*fArgs06(13) + iau06arr.as0i(i,14)*fArgs06(14);
        ssum4 = ssum4 + iau06arr.as0(i,1)*sin(tempval) + iau06arr.as0(i,2)*cos(tempval);
    end

    s = 0.000094 + 0.00380865 * ttt - 0.00012268 * ttt2  ...
        - 0.07257411 * ttt3 + 0.00002798 * ttt4 + 0.00001562 * ttt5; % ...
    %            + 0.00000171*ttt*sin(fArgs(5)) + 0.00000357*ttt*cos(2.0*fArgs(5)) ...
    %            + 0.00074353*ttt2*sin(fArgs(5)) + 0.00005691*ttt2*sin(2.0*(fArgs(3)-fArgs(4)+fArgs(5))) ...
    %            + 0.00000984*ttt2*sin(2.0*(fArgs(3)+fArgs(5))) - 0.00000885*ttt2*sin(2.0*fArgs(5));
    s = -x*y*0.5 + s*convrt + ssum0 + ssum1*ttt + ssum2*ttt2 + ssum3*ttt3 + ssum4*ttt4;  % rad

    if iauhelp == 'y'
        fprintf(1,'s %14.12f  %14.12f  %14.12f  %14.12f  %14.12f \n',ssum0/deg2rad,ssum1/deg2rad,ssum2/deg2rad,ssum3/deg2rad,ssum4/deg2rad );
    end

end