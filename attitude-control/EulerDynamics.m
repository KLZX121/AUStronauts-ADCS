function omegaDot = EulerDynamics(T, I)
%EulerDynamics Euler Rigid Body Dynamics Equation
%   omegaDot = EulerDynamics(T, omega0, I, dt)
%   
%   Inputs
%   T           (3, 1)  Total applied torque (N m)
%   I           (3, 3)  Moment of Inertia tensor (kg m^2)
%   
%   Outputs
%   omegaDot    (3, 1)  Angular velocity ODE (rad/s)
%   
%   TODO: account for time varying I (subsystem deployment) - update euler
%         equation
    
% Euler's equation (rigid body dynamics)
omegaDot = @(omega) I \ (T - cross(omega, I*omega));

end