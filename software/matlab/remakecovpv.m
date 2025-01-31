% ------------------------------------------------------------------------------
%
%                           function remakecovpv
%
%  takes propagated perturbed points from square root algorithm
%  and finds mean and covariance - position and velocity
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    pts            Array of propagated points from square root algorithm
%
%  outputs       :
%    cov             n_dim x n_dim covariance matrix (m)
%    yu             1 x n_dim mean vector (m)
%
%  locals        :
%    y              1 x n_dim mean shifted vector (m)
%    n_dim          dimension of vector
%    n_pts          total number of points
%
% ------------------------------------------------------------------------------

function [yu, cov] = remakecovpv(pts)

    % -------------------------  implementation   -----------------
    % initialize data & pre-allocate matrices
    [n_dim, n_pts] = size(pts);  % matrix is 12,6
    cov  = zeros(n_dim, n_dim);
    yu  = zeros(n_dim, 1);
    y   = zeros(n_dim, n_pts);
    tmp = zeros(n_dim, n_dim);

    % find mean
    for i = 1 : n_pts
        yu = yu + pts(1 : n_dim, i);
    end %  for i
    yu = yu / (n_pts);

    % find covariance
    for i = 1 : n_pts
        y(1:n_dim,i) = pts(1:n_dim,i) - yu(1:n_dim,1);
    end % for i
    cov = (y*y')/n_pts;

    %     cov = (tmp + tmp')/2; % ensures perfect symmetry
    %     cov1 = tmp;
    %     cov-cov1  % same


