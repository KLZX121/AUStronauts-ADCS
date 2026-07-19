function qDot = QKinematics(x, d)
%QKinematics Quaternion Kinematics ODE
%   qDot = QKinematics(x, d)
%   
%   Inputs
%   x       (:, 1)  State vector
%   d               Input parameters (struct)
%   .iQ     (1, 4)  Indices of quaternion in state vector
%   .iWSat  (1, 3)  Indices of sat ang vel in state vector
%
%   Outputs
%   qDot    (4, 1)  Quaternion derivative

% state variables
q = x(d.iQ);
w = x(d.iWSat);

Omega = @(w) [0 -w'; w -Skew(w)];

qDot = 0.5*Omega(w)*q;

end