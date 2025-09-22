
% ------------------------------------------------------------------------------
%
%                    function dcmtrvec
%
%  this function finds the vector resulting from a dcm transformation of another
%   vector.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                     range / units
%   dcm          - direction cosine matrix
%   vec          - input vector
%
%  outputs
%   vecout       - output vector
%
%  locals
%   none
%
%  coupling
%   none
%
%  references
%    xx
%
% [vecout]  = dcmtrvec ( dcm, vec ) ;
% ------------------------------------------------------------------------------

function [vecout] = dcmtrvec(dcm,vec)

       % ---- calculate output vector components
       vecout(1) = dcm(1,1)*vec(1) + dcm(1,2)*vec(2) + dcm(1,3)*vec(3);
       vecout(2) = dcm(2,1)*vec(1) + dcm(2,2)*vec(2) + dcm(2,3)*vec(3);
       vecout(3) = dcm(3,1)*vec(1) + dcm(3,2)*vec(2) + dcm(3,3)*vec(3);


