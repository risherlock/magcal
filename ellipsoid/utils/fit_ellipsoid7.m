function u = fit_ellipsoid7(x, y, z)
% Fits an ellipse with 7 parameters to the input sets of coordinates.
%
% Inputs:
%   x [nx1]: x-axis magnetometer measurements
%   y [nx1]: y-axis magnetometer measurements
%   z [nx1]: z-axis magnetometer measurements
%
% Output:
%   u [10x1]: [a,b,c,f,g,h,p,q,r,d], vector corresponding to the coefficients
%             of fitted ellipsoid where the general quadric surface given by
%             the equation:
%               ax2 + by2 + cz2 + 2fyz + 2gxz + 2hxy + 2px + 2qy + 2rz + d = 0
%             where, f = g = h = 0.
%
% risherlock (2025-02-10)

  % Design and covariance matrix
  n = length(x);
  D = [x.^2, y.^2, z.^2, 2*x, 2*y, 2*z, ones(n, 1)];
  S = D' * D;

  % Eigenvector to smallest eigenvalue
  [V, L] = eig(S);
  [~, idx] = min(diag(L));
  v = V(:, idx);

  % Extract ellipsoid parameters
  u = zeros(10, 1);
  u(1:3) = v(1:3); % a, b, c
  u(7:9) = v(4:6); % p, q, r
  u(10) = v(7); % d

  % Ensure positive definiteness of the ellipsoid matrix
  det_A = det([u(1), u(6), u(5); u(6), u(2), u(4); u(5), u(4), u(3)]);

  if det_A < 0
    u = -u;
  end
end
