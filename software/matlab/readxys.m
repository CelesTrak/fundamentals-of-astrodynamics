% -----------------------------------------------------------------------------
%
%                           function readXYS
%
%  this function initializes the xys iau2006 data. the input data files
%  are from processing the ascii files into a text file of xys calcualtion over
%  many years.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    xysarr      - array of xys data records                rad
%    xysLoc      - location for xys data file
%    infilename  - file name
%
%  outputs       :
%    xysarr      - array of xys data records                rad
%
%  locals        :
%
%  references    :
%
%   [xys06table] = readxys(infilename)
% -------------------------------------------------------------------------

function [xys06table] = readxys(infilename)

    % Define the full file path
    filepath = append(infilename, 'xysdata.dat');

    % Read the file directly into a table
    opts = delimitedTextImportOptions('NumVariables', 5);
    opts.DataLines = [1, Inf];
    opts.Delimiter = ' ';
    opts.VariableNames = {'jd', 'jdf', 'x', 'y', 's'};
    opts.VariableTypes = {'double', 'double', 'double', 'double', 'double'};

    xys06table = readtable(filepath, opts);

    % Compute derived columns directly in the table
    xys06table.mjd_tt = xys06table.jd + xys06table.jdf - 2400000.5;

end
