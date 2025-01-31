% ----------------------------------------------------------------------------
%
%                           function remakecovp
%
%  takes propagated perturbed points from square root algorithm
%  and finds mean and covariance - position only
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
%    references  :
%     alfano original code
%
% [yu, cov] = remakecovp(pts)
% ---------------------------------------------------------------------------

function [yu, cov] = remakecovp(pts)
    oo6 = 1.0 / 6.0;

    % -------------------------  implementation    ------------------
    % initialize data & pre-allocate matrices

    % find mean
    for i = 1 : 6
        for j = 1 : 3
            yu(j) = yu(j) + sigmapts(j, i);
        end
    end

    for j = 1 : 3
        yu(j) = yu(j) * oo6;
    end

    % find covariance
    for i = 1 : 6
        for j = 1 : 3
            y(j, i) = sigmapts(j, i) - yu(j);

            yt = mattransx(y, 3, 6);
            tmp = matmult(y, yt, 3, 6, 3);
            cov = matscale(tmp, 3, 3, oo6);
        end
    end

end % remakecovp
