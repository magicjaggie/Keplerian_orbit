function dy = ode_keplerian_orbit( ~, y, mu ) % y = []

r = y(1:3); % 3 param for position
v = y(4:6); % 3 param for velocity

R = sqrt(r(1).^2 + r(2).^2 + r(3).^2);

dy = [v(1);
    v(2);
    v(3);
    -mu/(R.^3)*r(1);
    -mu/(R.^3)*r(2);
    -mu/(R.^3)*r(3)];
end