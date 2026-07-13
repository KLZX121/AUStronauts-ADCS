function qDot = QKinematics(x, d)
%QKinematics Quaternion Kinematics ODE
%   qDot = QKinematics(q, omega)
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

Xi = @(q) [-q(2:4)'; q(1)*eye(3)+CrossProdM(q(2:4))];

qDot = 0.5*Xi(q)*w;

function M = CrossProdM(v)
    M = [0, -v(3), v(2);
         v(3), 0, -v(1);
         -v(2), v(1), 0;];
end

end