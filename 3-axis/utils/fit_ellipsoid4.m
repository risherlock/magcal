function u = fit_ellipsoid4(x, y, z)
% Fits an ellipse with 4 parameters to the input sets of coordinates.
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
%             where, f = g = h = p = q = r = 0.
%
% risherlock (2025-02-10)

  n = length(x);
  X = [2*x, 2*y, 2*z, ones(n, 1)];
  Y = x.^2 + y.^2 + z.^2;
  XTY = X' * Y;
  XTX = X' * X;
  v = -XTX \ XTY;

  % Extract ellipsoid parameters
  u = zeros(10, 1);
  u(1:3) = ones(size(u(1:3)));
  u(7:9) = v(1:3); % p, q, r
  u(10) = v(4); % d
end
