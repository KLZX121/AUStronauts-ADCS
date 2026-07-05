function qDot = QKinematics(q, omega)
%QKinematics Quaternion Kinematics ODE
%   qDot = QKinematics(q, omega)
%   
%   Inputs
%   q       (4, 1)  Quaternion
%   omega   (3, 1)  Angular velocity (rad/s)
%   
%   Outputs
%   qDot    (4, 1)  Quaternion derivative

Xi = @(q) [-q(2:4)'; q(1)*eye(3)+CrossProdM(q(2:4))];

qDot = 0.5*Xi(q)*omega;

function M = CrossProdM(v)
    M = [0, -v(3), v(2);
         v(3), 0, -v(1);
         -v(2), v(1), 0;];
end

end