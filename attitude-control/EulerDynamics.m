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
%   
%   TODO: account for time varying I (subsystem deployment) - update euler
%         equation

omegaDot = d.I \ (d.T - cross(omega, d.I*omega));

end