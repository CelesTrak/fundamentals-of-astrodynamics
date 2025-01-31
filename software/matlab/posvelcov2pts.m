% ------------------------------------------------------------------------------
%
%                           function posvelcov2pts
%
%  finds 12 sigma points given position and velocity and covariance information
%  using cholesky matrix square root algorithm, then propagates covariance points
%
%  author        : sal alfano      719-573-2600   31 mar 2011
%
%  revisions
%                  dave vallado   make single routine to simply find sigma
%                  points
%
%  inputs          description 
%    reci        - eci 3x1 position vector        km or m
%    veci        - eci 3x1 velocity vector        km or m
%    cov         - eci 6x6 covariance matrix      km or m
%
%  outputs       :
%    sigmapts    - structure of sigma points (6 x 12)
%
%  locals        :
%
%  sigmapts = posvelcov2pts(reci, veci, cov);
% ------------------------------------------------------------------------------

function sigmapts = posvelcov2pts(reci, veci, cov)

    % -------------------------  implementation   -----------------
    
    % initialize data & pre-allocate new points
    sigmapts = zeros(6, 12);  % array is 6 x 12

    % compute matrix square root of nP
    s = sqrt(6) * chol(cov,'lower');

    % perturb states, propagate to toff
    for i = 1:6  
        % ---- find positive direction vectors
        jj = (i-1)*2 + 1;  % incr this by 2
        sigmapts(1:3,jj) = reci + s(1:3,i); 
        sigmapts(4:6,jj) = veci + s(4:6,i);

        % ---- find negative direction vectors
        sigmapts(1:3,jj+1) = reci - s(1:3,i); 
        sigmapts(4:6,jj+1) = veci - s(4:6,i);
    end % for i

end

      

      
