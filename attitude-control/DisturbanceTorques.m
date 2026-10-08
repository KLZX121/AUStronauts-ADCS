function [T, surfData] = DisturbanceTorques(x, d, env, geom)
%DisturbanceTorques Calculates the four main external disturbance torques
%   [T, surfData] = DisturbanceTorques(x, d, env, geom)
%   
%   Aerodynamic and SRP forces are calculated using satellite geometry, as
%   defined by a collection of rectangular surfaces representing the
%   outside faces of the satellite. 
%
%   Inputs
%   x           (:, 1)      State vector
%   d           (struct)    Data struct
%   .iR         (1, 3)      Index of orbital position in state vector (km)
%   .iV         (1, 3)      Index of orbital velocity (km/s)
%   .iQ         (1, 4)      Index of attitude quaternion
%   .ISat       (3, 3)      Satellite MOI tensor (kg m^2)
%   env         (struct)    Environmental data struct
%   .rho        (1, 1)      Atmospheric density (kg m^-3)
%   .p          (1, 1)      Total Solar Radiation Pressure (N/m^2)
%   .s          (3, 1)      Unit Sun vector (ECI)
%   B           (3, 1)      Ambient magnetic field (T) (ECI)
%   geom        (struct)    Satellite geometry struct
%   .n_surfaces (1, 1)      Number of surfaces
%   .surfaces   (struct[])  Array of rectangular surface objects
%       .name       (str)   String label of surface
%       .normal     (3, 1)  Normal vector of surface
%       .area       (1, 1)  Area of surface (m^2)
%       .centroid   (3, 1)  Centroid location from centre of mass (m)
%       .Cd         (1, 1)  Coefficient of drag of the surface
%
%   Outputs
%   T       (struct)    Disturbance torque struct
%   .g      (3, 1)      Gravity gradient torque (N m)
%   .a      (3, 1)      Aerodynamic drag torque (N m)
%   .s      (3, 1)      SRP torque (N m)
%   .m      (3, 1)      Magnetic field torque (N m)
%   surfData   (struct)    Individual forces and torques for each surface
%   .a      (struct)    Aerodynamic drag
%       .f  (3, n)      Surface force (N)
%       .t  (3, n)      Surface torque (N m)
%   .s      (struct)    SRP drag
%       .f  (3, n)      Surface force (N)
%       .t  (3, n)      Surface torque (N m)

dcm = QToDCM(x(d.iQ));

%% gravity gradient
mu = 3.986e14; % m^3 s^-2

rb = dcm*x(d.iR); % km
rb = rb*1e3; % switch to m (delete once we update state vector to use m)

T.g = 3*mu/norm(rb)^5 * cross(rb, d.ISat*rb);

%% aero drag + SRP
% TODO: remove conversion to m
% TODO: add relative vel to atmosphere
v = x(d.iV)*1e3;

% convert velocity to body frame
v = dcm*v;
vMag = norm(v);
vHat = v ./ vMag;

% convert sun vec to body frame
s = dcm*env.s;
s = s ./ norm(s);

for i = 1:geom.n_surfaces
    surf = geom.surfaces(i);

    % projected area
    cosThetaAero = dot(surf.normal, vHat);
    ApAero =  heaviside(cosThetaAero) * surf.area * cosThetaAero;

    cosThetaSRP = dot(surf.normal, s);
    ApSRP = heaviside(cosThetaSRP) * surf.area * cosThetaSRP;

    % force
    fAero = (-vHat).*0.5*env.rho*vMag^2*ApAero*surf.Cd;

    fSRP = -env.p*ApSRP*((surf.Ca + surf.Crd)*s + (2/3*surf.Crd + 2*surf.Crs*cosThetaSRP)*surf.normal);


    % torque
    tAero = cross(surf.centroid, fAero);
    tSRP = cross(surf.centroid, fSRP);

    surfData.a.f(:, i) = fAero;
    surfData.a.t(:, i) = tAero;
    surfData.s.f(:, i) = fSRP;
    surfData.s.t(:, i) = tSRP;
end


T.a = sum(surfData.a.t, 2);

% TODO: don't apply in eclipse
T.s = sum(surfData.s.t, 2);


%% magnetic

T.m = [0; 0; 0];

end