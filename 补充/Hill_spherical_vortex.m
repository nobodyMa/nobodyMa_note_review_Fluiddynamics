% Hill's Spherical Vortex - Streamlines in meridional plane
% Parameters
R = 1.0;          % sphere radius
U = 1.0;          % free-stream velocity
A = 15*U/(2*R^2); % from Eq. (7.52b): U = (2/15)*A*R^2

% Grid
n = 200;
r = linspace(0, 3*R, n);
z = linspace(-3*R, 3*R, n);
[Rg, Zg] = meshgrid(r, z);   % Rg: radius, Zg: axial coordinate

% Compute stream function Psi
Psi = zeros(size(Rg));

% Inside sphere (r^2+z^2 < R^2)
inside = (Rg.^2 + Zg.^2) <= R^2;
Psi(inside) = -(A/10) * Rg(inside).^2 .* (R^2 - Zg(inside).^2 - Rg(inside).^2);

% Outside sphere (potential flow around a sphere)
outside = ~inside;
r_mag = sqrt(Rg(outside).^2 + Zg(outside).^2);
% Stream function for uniform flow + dipole: Psi = (U/2) * r^2 * (1 - R^3/r_mag^3)
% Here r is the cylindrical radius, r_mag is the spherical radius
Psi(outside) = (U/2) * Rg(outside).^2 .* (1 - R^3 ./ r_mag.^3);

% Plot contour lines (streamlines)
figure;
contour(Rg, Zg, Psi, 30, 'LineWidth', 1.2);
axis equal; grid on; hold on;

% Mark the sphere boundary
theta = linspace(0, 2*pi, 100);
x_circle = R * cos(theta);
y_circle = R * sin(theta);
plot(x_circle, y_circle, 'k--', 'LineWidth', 2);

% Labels and title
xlabel('r (cylindrical radius)');
ylabel('z (axial coordinate)');
title('Hill\'s Spherical Vortex: Streamlines in Meridional Plane');
legend('Streamlines (contours of \Psi)', 'Sphere boundary');
hold off;