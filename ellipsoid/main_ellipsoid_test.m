% Test of fit_ellipsoid10.m function
clc
clear 
close all

%%%%% Ellipsoid parameters for test %%%%%
% Semi principal axes
ax= 20; 
bx = 30; 
cx = 50;  

% Centre
xc = 10; 
yc = 100; 
zc = 10;  
centre = [xc,yc,zc];

% Rotation angles (3-2-1)
yaw = pi/4; 
pitch = pi/3; 
roll = pi/7;  

% Signal to noise ratio
SNR = 40;

%%%%% Dont change beyond this point %%%%%
% Generate test points
% x	=	a*cosu*sinv 	
% y	=	b*sinu*sinv	
% z	=	c*cosv
% u in [0,2pi) and v in [0,pi]
[u, v] = meshgrid(0:0.3:pi*2,0:0.3:pi);
x = ax*cos(u).*cos(v);
y = bx*cos(u).*sin(v);
z = cx*sin(u);

% Meshgrid to vector
x=x(:); y=y(:); z=z(:); 
xyz = [x y z];

% Rotate using DCM (321)
C = dcm321Euler(yaw,pitch,roll);
for i_iters = 1: length(x)
     new = C*xyz(i_iters,:)';
     xyz(i_iters,:) = new'; 
end

% Move centre after rotation
x = xc + xyz(:,1);
y = yc + xyz(:,2);
z = zc + xyz(:,3);

% Add noise to generated points
x = awgn(x,SNR,'measured');
y = awgn(y,SNR,'measured');
z = awgn(z,SNR,'measured');

% Ellipsoid fit algoritm
v = fit_ellipsoid10(x,y,z);

% Plot and camera settings
% set(gca,'NextPlot','add', 'Visible','off'); view(59,13); hold on;
plot3(x,y,z,'.','MarkerSize',15); hold on; % Input points
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Ellipsoid fitting');
grid on;
plot_ellipsoid(v);

function [C] = dcm321Euler(psi,theta,phi)
%%% DCM in terms of (3-2-1) Euler angles
% Maps a vector from inertial to body frame: V_b = C*V_i
%   
% Inputs:
%   psi     = The yaw rotation (radians)
%   theta   = The pitch rotation (radians)
%   phi     = The roll rotation (radians)
%
% Output:
%   C = Direction cosine matrix
%
% References:
%   Schaub - Analytical Mechanics of Space Systems (2nd ed.)(Pg.89)
%   Ardakani - Review of the 3-2-1 Euler Angles (2010)
%
% 2020/6/23

c1 = cos(psi);
s1 = sin(psi);
c2 = cos(theta);
s2 = sin(theta);
c3 = cos(phi);
s3 = sin(phi);

C = [   c2*c1,              c2*s1,              -s2;    ... 
        s3*s2*c1 - c3*s1,   s3*s2*s1 + c3*c1,   s3*c2;  ...    
        c3*s2*c1 + s3*s1,   c3*s2*s1 - s3*c1,   c3*c2  ];
end