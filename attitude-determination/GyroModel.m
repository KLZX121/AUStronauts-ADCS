function [omegaMeas, dOut] = GyroModel(omegaTrue, dt)
%GyroModel  Models a 3-axis rate gyro
%   [omegaMeas, d] = GyroModel(d, omegaTrue, dt)
%
%   Inputs 
%   omegaTrue   (3, 1)  True angular velocity (rad/s)
%   dt                  Timestep (s)
%
%   Outputs
%   omegaMeas   (3, 1)  Measured angular velocity (rad/s)
%   dOut                Output gyro data structure

% gyro data structure
dClean.measurement = 0;
dClean.scaleFactor = 1;
dClean.bias = 0;
dClean.temp = 0;
dClean.tempNominal = 0;
dClean.maxRate = 1e10;      % effectively no saturation
dClean.unitVector = [1;0;0];   % overwritten per axis below
dClean.resolution = 1e-12;     % effectively no quantisation
dClean.deadband = 0;
dClean.accelCoeff = [0 0];
dClean.biasTempCoeff = [0 0];
dClean.biasDayToDay  = 0;
dClean.biasRandomDrift = 0;
dClean.biasRandomWalk = 0;
dClean.scaleFactorDayToDay = 0;
dClean.scaleFactorTempCoeff = [0 0];
dClean.scaleFactorPositiveCoeff = [0 0];
dClean.scaleFactorNegativeCoeff = [0 0];
dClean.initialize = 1;

% set 3 orthogonal axes (x, y, z)
d3Gyros = repmat(dClean, 3, 1);
d3Gyros(1).unitVector = [1; 0; 0];
d3Gyros(2).unitVector = [0; 1; 0];
d3Gyros(3).unitVector = [0; 0; 1];

% model
dOut(1) = RateGyro(d3Gyros(1), omegaTrue, dt);
dOut(2) = RateGyro(d3Gyros(2), omegaTrue, dt);
dOut(3) = RateGyro(d3Gyros(3), omegaTrue, dt);

% extract measurement
omegaMeas = zeros(3,1);

for i = 1:3
    omegaMeas(i) = dOut(i).measurement * dOut(i).resolution;
end

end