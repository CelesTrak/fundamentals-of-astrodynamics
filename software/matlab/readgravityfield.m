% ----------------------------------------------------------------------------
%
%                           function readGravityField
%
%  this function reads and stores the gravity field coefficients. the routine is
%  can process either un-normalized or normalized coefficients. it's important to
%  set the normal parameter correctly as the gravData will contain c,s, or cNor, sNor
%  as needed. this should help minimize confusion with future use. because large
%  (>170) gravity fields cause normalization  problems, it's preferable to use
%  normalized gravity fields but you must also use a recursion that normalizes the
%  legendre functions. arrays are dimensioned from 0 because L,m go from 0. this may
%  cause difficulties in FOR, Matlab and languages that don't permit 0 array indices.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    fname       - filename for gravity field
%    normalized  - normalized or unnormalized                           'y', 'n'
%    startKtr    - first line of ata in the gravity file                0, 1, ...
%
%  outputs       :
%    order       - size of gravity field                                1..2160..
%    gravData    - structure containing the gravity field coefficients,
%                  radius of the Earth, and gravitational parameter
%
%  locals :
%    L, m        - degree and order indices
%    ktr         - starting line of data in each gravity file
%
%  coupling      :
%    none
%
%  references :
%    vallado       2022, 550-551
%
% [gravarr] = readgravityfield(fname, normal);
% ----------------------------------------------------------------------------

function [gravarr] = readgravityfield(fname, normal)

    % read the whole file at once into lines of an array
    filename = fname;
    [fid, message] = fopen(filename,'r');
    fclose(fid);
    %fprintf(1,'%s \n',filename);
    filedat =load(filename);

    for i = 1:size(filedat,1)
        L = filedat(i,1);
        m = filedat(i,2);
        if normal == 'y'
            gravarr.cNor(L+1, m+1) = filedat(i,3);
            gravarr.sNor(L+1, m+1) = filedat(i,4);
        else
            gravarr.c = filedat(i,3);
            gravarr.s = filedat(i,4);
        end
    end  % L

    % % Define the full file path
    % filepath =fname;
    %
    % % Read the file directly into a table
    % opts = delimitedTextImportOptions('NumVariables', 4);
    % opts.DataLines = [1, Inf];
    % opts.Delimiter = ' ';
    % opts.VariableNames = {'L', 'm', 'cNor', 'sNor'};
    % opts.VariableTypes = {'double', 'double', 'double', 'double'};
    %
    % gravarr = readtable(filepath, opts);


    % just read it in as it comes, doing the normalization or undoing
    % it gets messed up with large fields
    %sizeL = size(Larr);
    % i = 1;
    % while (Larr(i) < degsize)
    %     %for i = 1:degsize   %sizeL(1)
    %     % find normalized or unnormalized depending on which is already in file
    %     % note that above n = 170, the factorial will return 0, thus affecting the results!!!!
    %     L = Larr(i);
    %     m = marr(i);
    %     Lp1 = Larr(i) + 1;
    %     mp1 = marr(i) + 1;
    %
    %     if (m == 0)
    %         conv = sqrt( 2 * L + 1 );
    %         gravarr.conv(Lp1, mp1) = conv;
    %     else
    %         conv = sqrt( (factorial(L - m) * 2 * (2 * L + 1)) / factorial(L + m) );
    %         gravarr.conv(Lp1, mp1) = conv;
    %     end
    %     if normal == 'n'
    %         temp = 1.0 / conv;
    %         gravarr.cNor(Lp1, mp1) = temp * gravarr.c(i);
    %         gravarr.sNor(Lp1, mp1) = temp * gravarr.s(i);
    %         gravarr.c(Lp1, mp1) = gravarr.c(i);
    %         gravarr.s(Lp1, mp1) = gravarr.s(i);
    %     else
    %         gravarr.c(Lp1, mp1) = conv * gravarr.cNor(i);
    %         gravarr.s(Lp1, mp1) = conv * gravarr.sNor(i);
    %         gravarr.cNor(Lp1, mp1) = gravarr.cNor(i);
    %         gravarr.sNor(Lp1, mp1) = gravarr.sNor(i);
    %     end
    %     i = i + 1;
    % end  % while

end