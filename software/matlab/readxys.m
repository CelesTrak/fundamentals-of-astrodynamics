% -----------------------------------------------------------------------------
%
%                           function readxys
%
%  this function initializes the xys iau2006 iau data. the input data files
%  are from processing the ascii files into a text file of xys calcualtion over
%  many years. the jd and jdf for the start of the file are in the first
%  records. 
%  Note: Read using a table instead of a structure for efficiency. all data 
%  manipulations remain vectorized and efficient.
%
%  author        : david vallado           davallado@gmail.com   22 jan 2018
%
%  inputs          description                           range / units
%    infilename      - location for xys data file       
%
%  outputs       :
%    xys06table  - array of xys data records             rad
% 
%  locals        :
%
%  references    :
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
