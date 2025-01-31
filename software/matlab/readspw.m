% -----------------------------------------------------------------------------
%
%                           function readspw.m
%
%  this function reads the spw coefficients.
%
%  author        : david vallado                              29 nov 2024
%
%  revisions
%                -
%
%  inputs          description                           range / units
%    infilename  - filename
%
%  outputs       :
%    spwarr      - array of spw values
%
%  locals        :
%
%  coupling      :
%    none.
%
%  references    :
%    none.
%
% [spwarr, updDate] = readspw(spwFileName)
% -----------------------------------------------------------------------------

function [spwarr, updDate] = readspw(spwFileName)

    spwarr = struct('kparray',zeros(28000,8), 'sumkp',zeros(28000), ...
        'aparray',zeros(28000,8), 'avgap',zeros(28000), ...
        'adjf107',zeros(28000), 'adjctrf81',zeros(28000), 'adjlstf81',zeros(28000), ...
        'obsf107',zeros(28000), 'obsctrf81',zeros(28000), 'obslstf81',zeros(28000), ...
        'mjd',zeros(28000) );
    
    updDate = '';
    infile = fopen(spwFileName, 'r');

    while (~feof(infile))
        longstr = fgets(infile);
        while ( (longstr(1) == '#') && (feof(infile) == 0) )
            longstr = fgets(infile);
    
            if (contains(longstr, 'UPDATED'))
                updDate = longstr(8:end);
            end
    
            if (contains(longstr, 'NUM_OBSERVED_POINTS'))
                numrecsobs = str2double(longstr(20:end));
                longstr = fgets(infile);

                % ---- process observed records only
                % 1957 10 02 1700 20 37 37 17 17 27 23 17 30 203  22  22   6   6  12   9   6  15  12 0.7 3 234 253.5 6  265.7 231.4 253.3 267.4 231.7 233.3 233.0 215.9 211.6
                for ktr = 1: numrecsobs
                    longstr = fgets(infile);
                    [linearr] = sscanf( longstr, '%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %d %d %f %d %f %f %f %f %f');
                    
                    %spwarr(ktr).yr = linearr(1); 
                    %spwarr(ktr).mon = linearr(2);
                    %spwarr(ktr).day = linearr(3);
                    %spwarr(ktr).brsn = linearr(4);
                    %spwarr(ktr).nd = linearr(5);
                    spwarr(ktr).kparr(1) = linearr(6);
                    spwarr(ktr).kparr(2) = linearr(7);
                    spwarr(ktr).kparr(3) = linearr(8);
                    spwarr(ktr).kparr(4) = linearr(9);
                    spwarr(ktr).kparr(5) = linearr(10);
                    spwarr(ktr).kparr(6) = linearr(11);
                    spwarr(ktr).kparr(7) = linearr(12);
                    spwarr(ktr).kparr(8) = linearr(13);
                    spwarr(ktr).sumkp = linearr(14);
                    spwarr(ktr).aparr(1) = linearr(15);
                    spwarr(ktr).aparr(2) = linearr(16);
                    spwarr(ktr).aparr(3) = linearr(17);
                    spwarr(ktr).aparr(4) = linearr(18);
                    spwarr(ktr).aparr(5) = linearr(19);
                    spwarr(ktr).aparr(6) = linearr(20);
                    spwarr(ktr).aparr(7) = linearr(21);
                    spwarr(ktr).aparr(8) = linearr(22);
                    spwarr(ktr).avgap = linearr(23);
                    %spwarr(ktr).cp = linearr(24);
                    %spwarr(ktr).c9 = linearr(25);
                    %spwarr(ktr).isn = linearr(26);
                    spwarr(ktr).adjf107 = linearr(27);
                    %spwarr(ktr).q = linearr(28);
                    spwarr(ktr).adjctrf81 = linearr(29);
                    spwarr(ktr).adjlstf81 = linearr(30);
                    spwarr(ktr).obsf107 = linearr(31);
                    spwarr(ktr).obsctrf81 = linearr(32);
                    spwarr(ktr).obslstf81 = linearr(33);
                    %spwarr(ktr).adjf30 = linearr(34);
                    %spwarr(ktr).obsf30 = linearr(35);
                    %spwarr(ktr).adjlstf30f81 = linearr(36);
                    %spwarr(ktr).obslstf30f81 = linearr(37);
                    [jd, jdFrac] = jday(linearr(1), linearr(2), linearr(3), 0, 0, 0.0);
                    spwarr(ktr).mjd = jd + jdFrac - 2400000.5;
    
                end  % for through observed
                
            end
       
            % ---- process predicted records
            if (contains(longstr, 'NUM_DAILY_PREDICTED_POINTS'))
                numrecsobs = str2double(longstr(20:end));
                longstr = fgets(infile);
                % ---- process predicted records only
                % 1962 01 01 37665 -0.012700  0.213000  0.0326338  0.0017230  0.064261  0.006067  0.000000  0.000000   2 
                for i = 1: numrecsobs
                    longstr = fgets(infile);
                    [linearr] = sscanf( longstr, '%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %d %d %f %d %f %f %f %f %f');
                    
                    %spwarr(ktr+i).yr = linearr(1); 
                    %spwarr(ktr+i).mon = linearr(2);
                    %spwarr(ktr+i).day = linearr(3);
                    %spwarr(ktr+i).brsn = linearr(4);
                    %spwarr(ktr+i).nd = linearr(5);
                    spwarr(ktr+i).kparr(1) = linearr(6);
                    spwarr(ktr+i).kparr(2) = linearr(7);
                    spwarr(ktr+i).kparr(3) = linearr(8);
                    spwarr(ktr+i).kparr(4) = linearr(9);
                    spwarr(ktr+i).kparr(5) = linearr(10);
                    spwarr(ktr+i).kparr(6) = linearr(11);
                    spwarr(ktr+i).kparr(7) = linearr(12);
                    spwarr(ktr+i).kparr(8) = linearr(13);
                    spwarr(ktr+i).sumkp = linearr(14);
                    spwarr(ktr+i).aparr(1) = linearr(15);
                    spwarr(ktr+i).aparr(2) = linearr(16);
                    spwarr(ktr+i).aparr(3) = linearr(17);
                    spwarr(ktr+i).aparr(4) = linearr(18);
                    spwarr(ktr+i).aparr(5) = linearr(19);
                    spwarr(ktr+i).aparr(6) = linearr(20);
                    spwarr(ktr+i).aparr(7) = linearr(21);
                    spwarr(ktr+i).aparr(8) = linearr(22);
                    spwarr(ktr+i).avgap = linearr(23);
                    %spwarr(ktr+i).cp = linearr(24);
                    %spwarr(ktr+i).c9 = linearr(25);
                    %spwarr(ktr+i).isn = linearr(26);
                    spwarr(ktr+i).adjf107 = linearr(27);
                    %spwarr(ktr+i).q = linearr(28);
                    spwarr(ktr+i).adjctrf81 = linearr(29);
                    spwarr(ktr+i).adjlstf81 = linearr(30);
                    spwarr(ktr+i).obsf107 = linearr(31);
                    spwarr(ktr+i).obsctrf81 = linearr(32);
                    spwarr(ktr+i).obslstf81 = linearr(33);
                    %spwarr(ktr+i).adjf30 = linearr(34);
                    %spwarr(ktr+i).obsf30 = linearr(35);
                    %spwarr(ktr+i).adjlstf30f81 = linearr(36);
                    %spwarr(ktr+i).obslstf30f81 = linearr(37);
                    [jd, jdFrac] = jday(linearr(1), linearr(2), linearr(3), 0, 0, 0.0);
                    spwarr(ktr+i).mjd = jd + jdFrac - 2400000.5;
                end  % for through observed
               
            end
    
        end  % while through file

    end


%     opts = detectImportOptions('bigfile.txt'); 
% opts.DataLines = [3 8];
% opts.VariableNames = {'Timestamp','Temp',...
%                       'Humidity','Wind','Weather'};
% T_first = readtable('bigfile.txt',opts) 
% 
