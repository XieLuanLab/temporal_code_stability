function [x,y,z]=ellipsoid_Hanlin(meas,cov_matrix,Nsd,Npt)
% Find the eigenvalues and eigenvectors
[V,D] = eig(cov_matrix);

% Create a unit sphere
[phi, theta] = meshgrid(linspace(0, 2*pi, 2*Npt), linspace(0, pi, Npt));
x = cos(phi).*sin(theta);
y = sin(phi).*sin(theta);
z = cos(theta);

% Apply the eigenvectors and eigenvalues to the unit sphere
for i = 1:numel(x)
    vec = [x(i); y(i); z(i)];
    vec = V*Nsd*sqrt(D)*vec;
    x(i) = vec(1);
    y(i) = vec(2);
    z(i) = vec(3);
end
x=x+meas(1);
y=y+meas(2);
z=z+meas(3);