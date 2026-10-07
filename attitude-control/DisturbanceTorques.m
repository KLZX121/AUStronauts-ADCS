function [Tg, Ta, Ts, Tm] = DisturbanceTorques(x, d, s, B)
%DisturbanceTorques Calculates the four main external disturbance torques
%   [Tg, Ta, Ts, Tm] = DisturbanceTorques(x, d, s, B)
%   
%   Inputs
%   x       (:, 1)      State vector
%   d       (struct)    Data struct
%   .iR     (1, 3)      Index of orbital position in state vector (km)
%   .iV     (1, 3)      Index of orbital velocity (km/s)
%   .iQ     (1, 4)      Index of attitude quaternion
%   .ISat   (3, 3)      Satellite MOI tensor (kg m^2)
%   s       (3, 1)      Unit Sun vector (ECI)
%   B       (3, 1)      Ambient magnetic field (T) (ECI)
%
%   Outputs
%   Tg      (3, 1)      Gravity gradient torque (N m)
%   Ta      (3, 1)      Aerodynamic drag torque (N m)
%   Ts      (3, 1)      SRP torque (N m)
%   Tm      (3, 1)      Magnetic field torque (N m)

Att = QToDCM(x(d.iQ));


% gravity gradient
mu = 3.986e14; % m^3 s^-2

rb = Att*x(d.iR); % km
rb = rb*1e3; % switch to m (delete once we update state vector to use m)

Tg = 3*mu/norm(rb)^5 * cross(rb, d.ISat*rb);

% centre of pressure
rcp = ones(3, 1).*randn(3, 1);
rcp = (0.01/norm(rcp)) .*rcp; % make norm 0.01 m (1 cm)

% drag
Cd = 2;
rho = 1e-11; % kg m^-3
A = 0.1; % m^2

vb = x(d.iV)*1e3; % switch to km/s

Fa = (-vb./norm(vb))*0.5*rho*Cd*A*norm(vb)^2;
Ta = cross(rcp, Fa);

% srp (doesnt apply in eclipse)
Phi = 1367e3; % W m^-2
c = 3e8; % m s^-1
q = 0.6;
theta = abs(acos(dot(vb./norm(vb), s)));

s = Att*s;

Fs = -s.*(Phi/c * A * (1+q) * cos(theta));
Ts = cross(rcp, Fs);

% magnetic
D = randn(3, 1);
D = 0.1/norm(D) * D; % A m^2

B = Att*B;

Tm = cross(D, B);

end