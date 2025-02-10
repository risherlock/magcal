% Batch 3-axis compass calibration using least squares ellipsoid fitting
% 2020/06/03

clc
clear
close all

addpath('utils');

% Import raw magnetometer readings
raw = importdata('raw_data.csv');
x = raw(:,1);
y = raw(:,2);
z = raw(:,3);

% Ellipsoid fit calibration
v = ellipsoid_fit(x, y, z);
[Ainv, b, r, rmse] = magcal_ellipsoid(x, y, z);
m = [x, y, z]';
m_hat = Ainv * (m - b);

% Plot uncalibrated data
subplot(1,2,1);
scatter3(m_hat(1, :), m_hat(2, :), m_hat(3, :),'fill','MarkerFaceColor','red');
title({'Before magnetometer calibration','(Ellipsoid fitted)'});
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal;

% Plot calibrated data
subplot(1,2,2);
scatter3(m_hat(1, :), m_hat(2, :), m_hat(3, :),'fill','MarkerFaceColor','blue');
title({'After magnetometer calibration'});
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal;

% Print calibration matrices
fprintf('3D magnetometer calibration based on ellipsoid fitting');
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
fprintf('\nThe calibration equation:')
fprintf('\n      m_hat = Ainv * (m - b) \n\nwhere,')
fprintf('\n m = Three-axis magnetometer measurement [3x1]');
fprintf('\n m_hat = Calibrated magnetometer [3x1]');
fprintf('\n\nAinv =\n'); disp(Ainv);
fprintf('\nb =\n'); disp(b);
fprintf('\nr = '); disp(r);
fprintf('\nrmse = '); disp(rmse);
