%clc;
%clear all;


% Nb of revolutions:
k = 1;

% Orbit parameters
mu = 398600.44;
R_E = 6371; % Earth radius
a = 37000; % semi-major axis
e = 0.4;

% Time
T = 2*pi()*sqrt(a^3/mu); % period
tfin = T;
t0 = 0;
tspan = linspace(0,Tperiod, 100)'; % faster
%tspan = [t0,k*tfin]; % slower

% Initial conditions
p = a*(1-e^2);
r0 = a*(1-e); % = rp, radius at the perigee
v0 = sqrt(mu/p)*(1+e); % = vp, v at the perigee
Y0 = [r0 0 0 0 v0 0];

% Set options
options = odeset( 'RelTol', 1e-12, 'AbsTol', 1e-12 );

% Perform the integration
[ T, Y ] = ode113( @(t,y)ode_keplerian_orbit(t,y,mu), tspan, Y0, options);
r = [Y(:,1),Y(:,2),Y(:,3)];
v = [Y(:,4),Y(:,5),Y(:,6)];

% Angular momentum vector = h
h = cross(r,v);

% Eccentricity vector
ev = (cross(v,h))./mu-r/norm(r);

% Plot 
figure(1)
plot(T,r(:,1))
title('Trajectory')

% Draw planet
figure(2)
image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';

npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha = 1;
[x, y, z] = ellipsoid(0, 0, 0, R_E, R_E, R_E, npanels);
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
cdata = imread(image_file);
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');


% Aspect
hold on
axis equal
grid on
axis([-2*a,2*a,-2*a,2*a,-2*R_E,2*R_E])
xlabel('rx [L]');
ylabel('ry [L]');
zlabel('rz [L]');
title('Orbit');

% Earth aspect


% Plot the results
hold on
%quiver3(0,0,0,ev(1),ev(2),ev(3)) % display vectors
comet3(r(:,1),r(:,2),r(:,3))
n = 10; % n_quiver
quiver3(r(1:n:end,1),r(1:n:end,2),r(1:n:end,3),v(1:n:end,1),v(1:n:end,2),v(1:n:end,3)) % every n values of r or v


