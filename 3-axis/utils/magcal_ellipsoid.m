function [Ainv, b, r, srmse] = magcal_ellipsoid(x, y, z)
% Computes soft-iron correction matrix Ainv and hard-iron bias vector b.
%
% Inputs:
%   x [nx1]: x-axis magnetometer measurements
%   y [nx1]: y-axis magnetometer measurements
%   z [nx1]: z-axis magnetometer measurements
%
% Outputs:
%   Ainv   [3x3]: Soft-iron correction matrix
%   b      [3x1]: Hard-iron correction bias
%   r      [1x1]: Radius of corrected sphere (magnitude of magnetic field)
%   srmse [1x1]: Scaled Root Mean Squared Error between raw data and fitted ellipsoid
%                (scaled by 1/(2*r^2) to normalize relative to the sphere's radius)
%
% Calibration equation:
%   m_hat = Ainv * (m - b)
%
% risherlock (2025-02-09)

  % Get ellipsoid coefficients v = [a,b,c,f,g,h,p,q,r,d]
  v = ellipsoid_fit(x, y, z);

  % Ellipsoid in matrix form: Ax + k = 0
  A = [v(1), v(6), v(5);
       v(6), v(2), v(4);
       v(5), v(4), v(3)];
  k = [v(7); v(8); v(9)];
  det_A = det(A);

  % Soft-iron correction matrix
  Ainv = sqrtm(A ./ nthroot(det_A, 3));

  % Hard-iron bias vector
  b = -A \ k;

  % Magnitude of magnetic field intensity (normalized radius)
  sum_terms = A(1,1)*b(1)^2 + 2*A(1,2)*b(1)*b(2) + 2*A(1,3)*b(1)*b(3) + ...
              A(2,2)*b(2)^2 + 2*A(2,3)*b(2)*b(3) + A(3,3)*b(3)^2 - v(10);
  B = sqrt(abs(sum_terms));
  r = B / sqrt(nthroot(det_A, 3));

  %  Scaled RMSE: Deviation from a perfect sphere
  m = [x, y, z]; n = length(x);
  m_hat = (Ainv * (m.' - b)).';
  r_sq = sum(m_hat.^2, 2);
  err = r_sq - r^2;
  rmse = sqrt(mean(err.^2));
  srmse = real(rmse / (2 * r^2)); % Scaled RMSE
end
