function [] = plot_sphere(pc, r)
% Plots sphere
%
% Input:
%   pc = A point [xc, yc, zc] representating centre of sphere.
%   r = Radius of sphere.
%
% 2020/06/04

% Circle parametric equation
theta = linspace(0,2*pi,20);
phi = linspace(0,pi,20);
[theta,phi] = meshgrid(theta,phi);
rho = r;

x = pc(1) + rho*sin(phi).*cos(theta);
y = pc(2) + rho*sin(phi).*sin(theta);
z = pc(3) + rho*cos(phi);

m = mesh(x,y,z);
set(m,'FaceAlpha',0.4);
axis equal; hold on; grid on;
end
