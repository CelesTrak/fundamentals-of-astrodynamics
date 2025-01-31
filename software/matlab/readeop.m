% -----------------------------------------------------------------------------
%
%                           function readeop.m
%
%  this function reads the eop coefficients. the start mjd is contained in
%  mjd(1).
%
%  author        : david vallado                  719-573-2600   12 dec 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    infilename  - filename
%
%  outputs       :
%    eoparr      - array of eop values
%
%  locals        :
%
%  coupling      :
%    none.
%
%  references    :
%    none.
%
% [eoparr] = readeop(eopFileName)
% -----------------------------------------------------------------------------

function [eoparr] = readeop(eopFileName)

    eoparr = struct('xp',zeros(25000), ... % 'xperr',zeros(25000), ...
        'yp',zeros(25000), ...  % 'yperr',zeros(25000), ...
        'dut1',zeros(25000), ... % 'dut1err',zeros(25000), ...
        'lod',zeros(25000), ...
        'ddpsi',zeros(25000), 'ddeps',zeros(25000), ...
        'dx',zeros(25000), 'dy',zeros(25000), ...
        'mjd',zeros(25000), 'dat',zeros(25000) );

    infile = fopen(eopFileName, 'r');

    while ( (feof(infile) == 0) )
        longstr = fgets(infile);

        % if (contains(longstr, 'UPDATED'))
        %     eoparr.updDate = longstr(8:end);
        % end

        if (contains(longstr, 'NUM_OBSERVED_POINTS'))
            numrecsobs = str2double(longstr(20:end));
            longstr = fgets(infile);
            % ---- process observed records only
            % 1962 01 01 37665 -0.012700  0.213000  0.0326338  0.0017230  0.064261  0.006067  0.000000  0.000000   2
            for ktr = 1: numrecsobs
                longstr = fgets(infile);
                %eoparr(ktr).year = str2double(longstr(1:4));
                %eoparr(ktr).mon = str2double(longstr(5:8));
                %eoparr(ktr).day = str2double(longstr(8:11));
                eoparr(ktr).mjd = str2double(longstr(11:17));
                eoparr(ktr).xp = str2double(longstr(17:27));
                eoparr(ktr).xerr = 0.0;
                eoparr(ktr).yp = str2double(longstr(27:37));
                eoparr(ktr).yerr = 0.0;
                eoparr(ktr).dut1 = str2double(longstr(37:48));
                eoparr(ktr).dut1err = 0.0;
                eoparr(ktr).lod = str2double(longstr(48:59));
                eoparr(ktr).ddpsi = str2double(longstr(59:69));
                eoparr(ktr).ddeps = str2double(longstr(69:79));
                eoparr(ktr).dx = str2double(longstr(79:89));
                eoparr(ktr).dy = str2double(longstr(89:99));
                eoparr(ktr).dat = str2double(longstr(99:103));

            end  % for through observed

        end

        % ---- process predicted records
        if (contains(longstr, 'NUM_PREDICTED_POINTS'))
            numrecsobs = str2double(longstr(21:end));
            longstr = fgets(infile);
            % ---- process predicted records only
            % 1962 01 01 37665 -0.012700  0.213000  0.0326338  0.0017230  0.064261  0.006067  0.000000  0.000000   2
            for i = 1: numrecsobs
                longstr = fgets(infile);
                %eoparr(ktr+i).year = str2double(longstr(1:4));
                %eoparr(ktr+i).mon = str2double(longstr(5:8));
                %eoparr(ktr+i).day = str2double(longstr(8:11));
                eoparr(ktr+i).mjd = str2double(longstr(11:17));
                eoparr(ktr+i).xp = str2double(longstr(17:27));
                eoparr(ktr+i).xerr = 0.0;
                eoparr(ktr+i).yp = str2double(longstr(27:37));
                eoparr(ktr+i).yerr = 0.0;
                eoparr(ktr+i).dut1 = str2double(longstr(37:48));
                eoparr(ktr+i).dut1err = 0.0;
                eoparr(ktr+i).lod = str2double(longstr(48:59));
                eoparr(ktr+i).ddpsi = str2double(longstr(59:69));
                eoparr(ktr+i).ddeps = str2double(longstr(69:79));
                eoparr(ktr+i).dx = str2double(longstr(79:89));
                eoparr(ktr+i).dy = str2double(longstr(89:99));
                eoparr(ktr+i).dat = str2double(longstr(99:103));

            end  % for through predicted

        end

    end  % while through data portion of file


