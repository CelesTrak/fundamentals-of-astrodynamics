% -----------------------------------------------------------------------------
%
%                           function readjplde
%
%  this function initializes the jpl planetary ephemeris data. the input
%  data files are from processing the ascii files into a text file of sun
%  and moon positions.
%
%  author        : david vallado                  719-573-2600   22 jan 2018
%
%  revisions
%
%  inputs          description                    range / units
%
%
%
%
%  outputs       :
%    jpldearr    - array of jplde data records
%    jdjpldestart- julian date of the start of the jpldearr data
%
%  locals        :
%                -
%
%  coupling      :
%
%  references    :
%
%
%  [jpldearr] = readjplde(infilename)
%  -------------------------------------------------------------------------- */

function [jpldetable] = readjplde(infilename)

    % read the whole file at once into lines of an array
    filename = infilename;
    [fid, message] = fopen(filename,'r');
    fclose(fid);
    %fprintf(1,'%s \n',filename);
    filedat =load(filename);

    %load data into x y z arrays
    year = filedat(:,1);
    mon = filedat(:,2);
    day = filedat(:,3);
    hr = filedat(:,4);
    jpldetable.rsun1 = filedat(:,5);
    jpldetable.rsun2 = filedat(:,6);
    jpldetable.rsun3 = filedat(:,7);
    jpldetable.rsmag = filedat(:,8);
    jpldetable.rmoon1 = filedat(:,10);
    jpldetable.rmoon2 = filedat(:,11);
    jpldetable.rmoon3 = filedat(:,12);

    [jdtdb, jdtdbf] = jday(year, mon, day, hr, 0, 0.0);
    jpldetable.mjd = jdtdb + jdtdbf - 2400000.5;


    % % Define the full file path
    % filepath = infilename;
    % 
    % % Read the file directly into a table
    % opts = delimitedTextImportOptions('NumVariables', 12);
    % opts.DataLines = [1, Inf];
    % opts.Delimiter = {' ','\t','\b'};
    % %opts.Whitespace  = ' ';
    % opts.ConsecutiveDelimitersRule = 'join';
    % opts.VariableNames = {'year', 'mon', 'day', 'hr', 'rsun1', 'rsun2', 'rsun3', 'rsmag', 'rmoon1', 'rmoon2', 'rmoon3', 'rmmag'};
    % opts.VariableTypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
    % 
    % jpldetable = readtable(filepath, opts);
    % 
    % % Compute derived columns directly in the table 
    % [jdtdb, jdtdbf] = jday(jpldetable.year, jpldetable.mon, jpldetable.day, jpldetable.hr, 0, 0.0);
    % jpldetable.mjd = jdtdb + jdtdbf - 2400000.5;
    
end  %  initjplde
