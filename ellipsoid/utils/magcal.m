function [Ainv, b, r, err] = magcal(x, y, z)
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
%   srmse  [1x1]: Scaled Root Mean Squared Error between raw data and fitted ellipsoid
%                (scaled by 1/(2*r^2) to normalize relative to the sphere's radius)
%
% Calibration equation:
%   m_hat = Ainv * (m - b)
%
% risherlock (2025-02-09)

  [A, Ainv, b, r, err, ispd] = get_all_params(x, y, z, 4);
  [A7, Ainv7, b7, r7, err7, ispd7] = get_all_params(x, y, z, 7);
  isrealA7 = isreal(A7) || all(imag(A7) == 0(:));

  if ispd7 && isrealA7 && (err7 < err)
    Ainv = real(Ainv7); b = b7; r = r7; err = err7;
  end

  [A10, Ainv10, b10, r10, err10, ispd10] = get_all_params(x, y, z, 10);
  isrealA10 = isreal(A10) || all(imag(A10) == 0(:));

  if ispd10 && isrealA10 && (err10 < err)
    Ainv = real(Ainv10); b = b10; r = r10; err = err10;
  end
end

% Get calibration params, fit quality, and field intensity
function [A, Ainv, b, r, err, ispd] = get_all_params(x, y, z, idx)
  if idx == 4
    v = fit_ellipsoid4(x, y, z);
  elseif idx == 7
    v = fit_ellipsoid7(x, y, z);
  else
    v = fit_ellipsoid10(x, y, z);
  end

  [A, Ainv, b] = get_calibration_params(v);
  r = get_field_intensity(v, A, b);
  [err, ispd] = get_fit_quality(x, y, z, A, Ainv, b, r);
end

% Ellipsoid coefficients to calibration params
function [A, Ainv, b] = get_calibration_params(v)
  % Ellipsoid in matrix form: Ax + k = 0
  A = [v(1), v(6), v(5);
       v(6), v(2), v(4);
       v(5), v(4), v(3)];
  k = [v(7); v(8); v(9)];

  % Soft-iron correction matrix
  Ainv = sqrtm(A ./ nthroot(det(A), 3));

  % Hard-iron bias vector
  b = -A \ k;
end

% Magnitude of magnetic field intensity (normalized radius)
function [r] = get_field_intensity(v, A, b)
  terms =  A(1,1) * b(1)^2 + 2 * A(1,2) * b(1) * b(2)
         + A(2,2) * b(2)^2 + 2 * A(1,3) * b(1) * b(3)
         + A(3,3) * b(3)^2 + 2 * A(2,3) * b(2) * b(3) - v(10);
  r = sqrt(abs(terms)) / sqrt(nthroot(det(A), 3));
end

% Quality of ellipsoid fit
function [srmse, ispd] = get_fit_quality(x, y, z, A, Ainv, b, r)
  %  Scaled RMSE: Deviation from a perfect sphere
  m = [x, y, z]; n = length(x);
  m_hat = (Ainv * (m.' - b)).';
  r_sq = sum(m_hat.^2, 2);
  err = r_sq - r^2;
  rmse = sqrt(mean(err.^2));
  srmse = real(rmse / (2 * r^2));

  % Is positive definite?
  [~,p] = chol(A);
  ispd = (p == 0);
end
