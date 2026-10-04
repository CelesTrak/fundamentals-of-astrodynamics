% ------------------------------------------------------------------------------
%
%                            function decodealpha5
%
%  this function decodes the 5-character satellite number field of a tle. the
%    field is either the original five digits, or an alpha-5 field - a letter
%    giving the leading 10 to 33, followed by four digits. the letters i and o
%    are never used in an alpha-5 field as they are too easily read as 1 and 0,
%    so a field beginning with either, or one whose last four characters are not
%    digits, is not a satellite number and is refused.
%
%  inputs          description                    range / units
%    satnumStr   - satellite number field         '00005', 'E8493'
%
%  outputs       :
%    satnum      - satellite number, -1 if the field is not valid
%    valid       - true if the field was decoded, false if it is not valid
%
%  locals        :
%    alpha5      - the alpha-5 blocks, a to z with i and o left out
%    ch          - satellite number field, trimmed
%    index       - leading block of the alpha-5 field, 10 to 33
%
%  coupling      :
%    none.
%
%  references    :
%    space-track.org "what is alpha-5?"
%
% [satnum, valid] = decodealpha5( satnumStr );
% ----------------------------------------------------------------------------- }

function [satnum, valid] = decodealpha5 ( satnumStr )

    alpha5 = [10 11 12 13 14 15 16 17 0 18 19 20 21 22 0 23 24 25 26 27 28 29 30 31 32 33];

    satnum = -1;
    valid = false;

    ch = strtrim(char(satnumStr));
    if (isempty(ch) || (length(ch) > 5))
        return;
    end

    % ---- the original satellite number, all digits ----
    if (all(ch >= '0' & ch <= '9'))
        satnum = str2double(ch);
        valid = true;
        return;
    end

    % ---- an alpha-5 field, one letter and four digits ----
    if (length(ch) ~= 5)
        return;
    end
    if (~all(ch(2:5) >= '0' & ch(2:5) <= '9'))
        return;
    end
    first = upper(ch(1));
    if ((first < 'A') || (first > 'Z'))
        return;
    end

    index = alpha5(double(first) - double('A') + 1);
    % i and o are not part of the scheme
    if (index == 0)
        return;
    end

    satnum = index * 10000 + str2double(ch(2:5));
    valid = true;
