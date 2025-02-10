function [] = plot_ellipsoid(v)
% Plots ellipsoid from coefficients of general quadrics equation.
% Result of whole night discussion with Manjul (2020/6/8-9).
%
% Input:
%   v = [a,b,c,f,g,h,p,q,r,d], a vector corresponding to coefficients of
%       the general quadric surface given by equation,
%       ax2 + by2 + cz2 + 2fyz + 2gxz + 2hxy + 2px + 2qy + 2rz + d = 0.
%
% Notes:
% 1.General ellipsoid in matrix form: X'*M*X + [2p,2q,2r]*X + d = 0
%   where,
%       X = [x,y,z]
%       M = [a h g; h b f; g f c]
%
% 2.Diagonalizing above matrix M to M_ rotates the general ellipsoid such
%   that its axes aligns with new coordinate axes defined by eigenvectors.
%   New ellipsoid: X'*Q'*M*Q*X + [2p,2q,2r]*Q*X + d = 0
%   where,
%       Q  = Eigenvector matrix of M (DCM from old to new coordinates)
%       M_ = Q'*M*Q= New ellipsoid diagonal matrix
%
%
% 3.Semi principal axes (ax,by,cz) for ellipsoid with no rotation
%       ax = sqrt(p^2/a^2 + q^2/(a*b) + r^2/(a*c) - d/a);
%       by = sqrt(p^2/(a*b) + q^2/b^2 + r^2/(b*c) - d/b);
%       cz = sqrt(p^2/(a*c) + q^2/(b*c) + r^2/c^2 - d/c);
%   here, f = 0, g = 0 and h = 0
%
% Steps:
%   1) Diagonalize M matrix using eigenvector of M, M_ = Q'*M*Q
%   2) Find p,q,r of new matrix, [p_,q_,r_] = [p,q,r]*Q
%   3) Find semi principal axes ax_,ay_,az_ for f = g = h = 0
%   4) Generate points of ellipsoid using parametric equation
%   5) Rotate the points using DCM, Q
%   6) Plot the points
%
% 2020/6/8

% Unpack ellipsoid coefficients
a = v(1); b = v(2); c = v(3);
f = v(4); g = v(5); h = v(6);
p = v(7); q = v(8); r = v(9);
d = v(10);

% Coordinate frame transformation i.e diagonalize M
M =[a h g; h b f; g f c]; % Original ellipsoid matrix
[evec,~]=eig(M); % Compute eigenvectors matrix
Q = evec; % DCM = eigenvectors matrix
M_ = Q'*M*Q; % Diagonalize M

% Coefficients of the ellipsoid in new frame
% Note the ellipsoid is not rotating in this new frame so f, g and h = 0
pqr_ = [p,q,r]*Q;
a_ = M_(1,1);
b_ = M_(2,2);
c_ = M_(3,3);
p_ = pqr_(1);
q_ = pqr_(2);
r_ = pqr_(3);
d_ = d;

% Semi principal axes (Still no rotation)
ax_ = sqrt(p_^2/a_^2 + q_^2/(a_*b_) + r_^2/(a_*c_) - d_/a_);
bx_ = sqrt(p_^2/(a_*b_) + q_^2/b_^2 + r_^2/(b_*c_) - d_/b_);
cx_ = sqrt(p_^2/(a_*c_) + q_^2/(b_*c_) + r_^2/c_^2 - d_/c_);

% Centre of the ellipsoid
centre = M\[-p, -q, -r]';

% Generate ellipsoid points
theta = linspace(0,2*pi,20);
phi = linspace(0,pi,20);
[theta,phi] = meshgrid(theta,phi);
x_ = ax_*cos(theta).*cos(phi);
y_ = bx_*cos(theta).*sin(phi);
z_ = cx_*sin(theta);

%%% Note
% Till now we have extracted all the information from the
% input coefficients of ellipsoid. We generated the points of ellipsoid
% too. It must be remembered that these points are described in the
% coordinate represented by eigenvectors where there is no rotation of
% ellipsoid. Now all that remains is to rotate each points using DCM (Q)
% that we found above.

% Rotate ellipsoid
grid_size = length(x_);
x_ = x_(:); y_ = y_(:); z_ = z_(:);
for i_iters = 1:length(x_)
    after_rot = Q*[x_(i_iters), y_(i_iters), z_(i_iters)]';
    x_(i_iters) = centre(1) + after_rot(1);
    y_(i_iters) = centre(2) + after_rot(2);
    z_(i_iters) = centre(3) + after_rot(3);
end

% Vector to mesh
x_ = reshape(x_,grid_size,grid_size);
y_ = reshape(y_,grid_size,grid_size);
z_ = reshape(z_,grid_size,grid_size);

% Plot the ellipsoid corresponding to input coefficients
m = mesh(x_,y_,z_,'MarkerFaceColor','red', 'FaceAlpha', 0.4);
set(m,'facealpha',0.4);
hold on; axis equal;
end
