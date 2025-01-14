% ----------------------------------------------------------------------------
%
%                           function iau06in
%
%  this function initializes the matricies needed for iau 2006 reduction
%    calculations. the routine uses the files listed as inputs, but they are
%    are not input to the routine as they are static files.
%
%  author        : david vallado                  719-573-2600   16 jul 2004
%
%  revisions
%    dav 14 apr 11  update for iau2006 conventions
%
%  inputs          description                    range / units
%    none
%    iau00x.dat  - file for x coefficient
%    iau00y.dat  - file for y coefficient
%    iau00s.dat  - file for s coefficient
%    iau00n.dat  - file for nutation coefficients
%    iau00pl.dat notused - file for planetary nutation coefficients
%    iau00gs.dat - file for gmst coefficients
%
%  outputs       :
%    ax0        - real coefficients for x        rad
%    a0xi        - integer coefficients for x
%    ay0        - real coefficients for y        rad
%    a0yi        - integer coefficients for y
%    as0        - real coefficients for s        rad
%    a0si        - integer coefficients for s
%    apn         - real coefficients for nutation rad
%    apni        - integer coefficients for nutation
%    ape         - real coefficients for obliquity rad
%    apei        - integer coefficients for obliquity
%    agst        - real coefficients for gst      rad
%    agsti       - integer coefficients for gst
%
%  locals        :
%    convrt      - conversion factor to radians
%    i           - index
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado     2004, pg 205-219, 910-912
%
% [iau06arr] = iau06in(infilename)
% ----------------------------------------------------------------------------- }

function [iau06arr] = iau06in(infilename)

    % ------------------------  implementation   -------------------
    % " to rad
    convrtu= (0.000001*pi) /(180.0*3600.0);  % if micro arcsecond
    convrtm= (0.001*pi) /(180.0*3600.0);     % if milli arcsecond

    % ------------------------------
    %  note that since all these coefficients have only a single
    %  decimal place, one could store them as integres, and then simply
    %  divide by one additional power of ten. it would make memeory
    %  storage much smaller and potentially faster.
    % ------------------------------

    iau06arr = struct('ax0',zeros(1600, 2), 'ax0i',zeros(1600, 14), ...
        'ay0',zeros(1275, 2), 'ay0i',zeros(1275, 14), ...
        'as0',zeros(66, 2), 'as0i',zeros(66, 14), ...
        'agst',zeros(34, 2), 'agsti',zeros(34), ...
        'apn0',zeros(1358, 2), 'apn0i',zeros(1358, 14), ...
        'apl0',zeros(687, 2), 'apl0i',zeros(687, 14), ...
        'aapn0',zeros(678, 6), 'aapn0i',zeros(678, 5), ...
        'aapl0',zeros(687, 5), 'aapl0i',zeros(687, 14) );

    % xys values
    fid = fopen(append(infilename,'iau06xtab5.2.a.dat'));  % 1600
    %      1    -6844318.44        1328.67    0    0    0    0    1    0    0    0    0    0    0    0    0    0
    Dataarr = textscan( fid, '%d %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d ', 'Headerlines', 2, 'MultipleDelimsAsOne', true );
    fclose(fid);
    iau06arr.ax0(:,1)   = cell2mat(Dataarr(:,2));
    iau06arr.ax0(:,2)   = cell2mat(Dataarr(:,3));
    iau06arr.ax0i(:,1:14)   = cell2mat(Dataarr(:,4:17));
    for i=1:size(iau06arr.ax0)
        iau06arr.ax0(i,1)= iau06arr.ax0(i,1) * convrtu;  % rad
        iau06arr.ax0(i,2)= iau06arr.ax0(i,2) * convrtu;  % rad
    end


    fid = fopen(append(infilename,'iau06ytab5.2.b.dat')); % 1275
    %       1       1538.18     9205236.26    0    0    0    0    1    0    0    0    0    0    0    0    0    0
    Dataarr = textscan( fid, '%d %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d ', 'Headerlines', 2, 'MultipleDelimsAsOne', true );
    fclose(fid);
    iau06arr.ay0(:,1)   = cell2mat(Dataarr(:,2));
    iau06arr.ay0(:,2)   = cell2mat(Dataarr(:,3));
    iau06arr.ay0i(:,1:14)   = cell2mat(Dataarr(:,4:17));
    for i=1:size(iau06arr.ay0)
        iau06arr.ay0(i,1)= iau06arr.ay0(i,1) * convrtu;
        iau06arr.ay0(i,2)= iau06arr.ay0(i,2) * convrtu;
    end


    fid = fopen(append(infilename,'iau06stab5.2.d.dat'));  % 66
    %      1    -6844318.44        1328.67    0    0    0    0    1    0    0    0    0    0    0    0    0    0
    Dataarr = textscan( fid, '%d %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d ', 'Headerlines', 2, 'MultipleDelimsAsOne', true );
    fclose(fid);
    iau06arr.as0(:,1)   = cell2mat(Dataarr(:,2));
    iau06arr.as0(:,2)   = cell2mat(Dataarr(:,3));
    iau06arr.as0i(:,1:14)   = cell2mat(Dataarr(:,4:17));
    for i=1:size(iau06arr.as0)
        iau06arr.as0(i,1)= iau06arr.as0(i,1) * convrtu;
        iau06arr.as0(i,2)= iau06arr.as0(i,2) * convrtu;
    end


    fid = fopen(append(infilename,'iau06ansofa.dat'));  % 678
    %  0     0     0     0     1               -172064161.0  -174666.0   33386.0  92052331.0   9086.0  15377.0
    Dataarr = textscan( fid, '%d %d %d %d %d %f %f %f %f %f %f ', 'Headerlines', 2, 'MultipleDelimsAsOne', true );
    fclose(fid);
    iau06arr.aapn0(:,1:6)   = cell2mat(Dataarr(:,6:11));
    iau06arr.aapn0i(:,1:5)   = cell2mat(Dataarr(:,1:5));
    %     2000a nutation values old approach iau2000a
    for i=1:size(iau06arr.aapn0)
        iau06arr.aapn0(i,1)= iau06arr.aapn0(i,1) * convrtm;
        iau06arr.aapn0(i,2)= iau06arr.aapn0(i,2) * convrtm;
        iau06arr.aapn0(i,3)= iau06arr.aapn0(i,3) * convrtm;
        iau06arr.aapn0(i,4)= iau06arr.aapn0(i,4) * convrtm;
        iau06arr.aapn0(i,5)= iau06arr.aapn0(i,5) * convrtm;
        iau06arr.aapn0(i,6)= iau06arr.aapn0(i,6) * convrtm;
    end

    % nutation values old approach iau2003
    %filein = load('iau03n.dat');
    % 2000a planetary nutation values
    fid = fopen(append(infilename,'iau06anplsofa.dat'));  % 687
    %       0   0   0   0   0   0   0   8 -16   4   5   0   0   0                   1440           0           0           0
    Dataarr = textscan( fid, '%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d ', 'Headerlines', 2, 'MultipleDelimsAsOne', true );
    fclose(fid);
    iau06arr.apl0(:,1:4)   = cell2mat(Dataarr(:,15:18));
    iau06arr.apl0i(:,1:14)   = cell2mat(Dataarr(:,1:14));
    for i=1:size(iau06arr.apl0)
        iau06arr.apl0(i,1)= iau06arr.apl0(i,1) * convrtm;
        iau06arr.apl0(i,2)= iau06arr.apl0(i,2) * convrtm;
        iau06arr.apl0(i,3)= iau06arr.apl0(i,3) * convrtm;
        iau06arr.apl0(i,4)= iau06arr.apl0(i,4) * convrtm;
    end

    % tab5.3a.txt in IERS
    % nutation values planetary now included new iau2006
    %     1   -17206424.18        3338.60    0    0    0    0    1    0    0    0    0    0    0    0    0    0
    % nutation in longitude
    fid = fopen(append(infilename,'iau06nlontab5.3.a.dat'));  % 1358
    Dataarr = textscan( fid, '%d %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d ', 'Headerlines', 2, 'MultipleDelimsAsOne', true );
    fclose(fid);
    iau06arr.apn0(:,1:2)   = cell2mat(Dataarr(:,2:3));
    iau06arr.apn0i(:,1:14)   = cell2mat(Dataarr(:,4:17));
    for i=1:size(iau06arr.apn0)
        iau06arr.apn0(i,1)= iau06arr.apn0(i,1) * convrtm;
        iau06arr.apn0(i,2)= iau06arr.apn0(i,2) * convrtm;
    end


    % tab5.3b.txt in IERS
    % nutation in obliquity
    %        fileData = File.ReadAllLines(nutLoc + "iau06nobltab5.3.b.dat");  // 1056


    % nutation values planetary now included new iau2006
    %       load iau00n.dat;  % luni-solar
    %       apn  = iau00n(:,2:3);
    %       apni   = iau00n(:,4:17);
    %       for i=1:size(apn)
    %           apn(i,1)= apn(i,1) * convrtu;
    %           apn(i,2)= apn(i,2) * convrtu;
    %       end;
    %
    %       load iau00e.dat;  % planetary
    %       ape  = iau00n(:,2:3);
    %       apei   = iau00n(:,4:17);
    %       for i=1:size(ape)
    %           ape(i,1)= ape(i,1) * convrtu;
    %           ape(i,2)= ape(i,2) * convrtu;
    %       end;

    % tab5.2e.txt in IERS
    % gmst values
    % note - these are very similar to the first 34 elements of iau00s.dat,
    % but they are not the same.
    fid = fopen(append(infilename,'iau06gsttab5.2.e.dat'));
    Dataarr = textscan( fid, '%d %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d ', 'Headerlines', 2, 'MultipleDelimsAsOne', true );
    fclose(fid);
    iau06arr.agst(:,1:2)   = cell2mat(Dataarr(:,2:3));
    iau06arr.agsti(:,1:14)   = cell2mat(Dataarr(:,4:17));
    for i=1:size(iau06arr.agst)
        iau06arr.agst(i,1)= iau06arr.agst(i,1) * convrtu;
        iau06arr.agst(i,2)= iau06arr.agst(i,2) * convrtu;
    end


