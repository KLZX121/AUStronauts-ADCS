function omegaDot = EulerDynamics(omega, d, w)
%EulerDynamics Euler Rigid Body Dynamics Equation
%   omegaDot = EulerDynamics(omega, d, w)
%   
%   INPUT
%   omega           (3, 1)  Angular velocity of reaction wheels (rad/s) 
%   d               Input parameters (struct)
%           .T      (3, 1)  Total external applied torque (N m)
%           .TRWA   (3, 1)  Total reaction wheel torque (N m)
%           .I      (3, 3)  Moment of Inertia tensor of satellite (kg m^2)
%           .IRWA   (3, 3)  Moment of Inertia of Reaction Wheels (kg m^2)
%
%   w               (3, 1)  Angular velocity of satellite (rad/s)
%   
%   OUTPUT
%   omegaDot        (3, 1)  Angular velocity ODE (rad/s)
%   
%   TODO: account for time varying I (subsystem deployment) - update euler
%         equation

c    = omega; % current reaction wheel rates (rad/s^2)
hT   = d.I*w + d.IRWA*c; % find change in momentum (kg m^2/s^2)
omegaDot = d.I\(d.T - d.TRWA - cross(w, hT)); % satellite angular velocity (rad/s^2)

end




%{
old function

function omegaDot = EulerDynamics(omega, d)
%EulerDynamics Euler Rigid Body Dynamics Equation
%   omegaDot = EulerDynamics(omega, d)
%   
%   Inputs
%   omega   (3, 1)  Angular velocity (rad/s)
%   d               Input parameters (struct)
%           .T      (3, 1)  Total applied torque (N m)
%           .I      (3, 3)  Moment of Inertia tensor (kg m^2)
%   
%   Outputs
%   omegaDot    (3, 1)  Angular velocity ODE (rad/s)

omegaDot = d.I \ (d.T - cross(omega, d.I*omega));

end

%}