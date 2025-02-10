function [u] = ellipsoid_fit(x,y,z)
% Ellipsoid fitting algorithm
%
% Inputs:
%   x [nx1]: x coordinates of input data
%   y [nx1]: y coordinates of input data
%   z [nx1]: z coordinates of input data
%
% Output:
%   u [10x1]: [a,b,c,f,g,h,p,q,r,d], vector corresponding to the coefficients
%             of fitted ellipsoid where the general quadric surface given by
%             the equation:
%               ax2 + by2 + cz2 + 2fyz + 2gxz + 2hxy + 2px + 2qy + 2rz + d = 0.
%
% Source:
%   [1] Li - Least Square Ellipsoid Fitting (2004)
%
% 2020/06/03

  % Design matrix
  D = [x.^2, y.^2, z.^2, 2*y.*z, 2*x.*z, 2*x.*y, 2*x, 2*y, 2*z, ones(size(x))]';

  % Inverse of constraint matrix with k = 4, Eqn(7)
  invC = [0,   0.5, 0.5, 0,    0,    0;
          0.5, 0,   0.5, 0,    0,    0;
          0.5, 0.5, 0,   0,    0,    0;
          0,   0,   0,  -0.25, 0,    0;
          0,   0,   0,   0,   -0.25, 0;
          0,   0,   0,   0,    0,   -0.25];

  % Eqn(11)
  S = D * D';
  S11 = S(1:6, 1:6);   % 6X6
  S12 = S(1:6, 7:10);  % 6X4
  S22 = S(7:10, 7:10); % 4X4
  Sx = pinv(S22) * S12';

  % Eqn(14) and Eqn(15)
  M = invC * (S11 - S12 * Sx);
  [V, L] = eig(M);

  % Index of the 'only' positive eigenvalue.
  [max_ridx, ~] = max(L);
  [~, max_cidx] = max(max_ridx);

  % Eigenvector corresponding to max_cidx
  u1 = V(:, max_cidx);
  u2 = -Sx * u1;
  u = [u1', u2']';

  % Ellipsoid in matrix form: Ax + k = 0
  det_A = det([u(1), u(6), u(5); u(6), u(2), u(4); u(5), u(4), u(3)]);

  % Ensure positive definiteness
  if det_A < 0
    u = -u;
  end
end
