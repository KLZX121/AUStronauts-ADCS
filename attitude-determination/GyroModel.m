function wGyro = GyroModel(wTrue, gyro, dt)
%GyroModel  Models a 3-axis rate gyro
%   wGyro = GyroModel(wTrue, gyro, dt)
%
%   Inputs 
%   wTrue   (3, 1)      True angular velocity (rad/s)
%   gyro    (struct)    Gyro calibration parameters
%   .bias   (3, 1)      Gyro bias (rad/s)
%   .s      (3, 1)      Scale factors
%   .kU     (3, 1)      Upper misalignments
%   .kL     (3, 1)      Lower misalignments
%   .S      (3, 3)      Scale factor/misalignment matrix
%   .sigG   (1, 1)      1-sigma measurement noise
%   .sigB   (1, 1)      1-sigma bias walk
%   dt                  Timestep (s)
%
%   Outputs
%   wGyro   (3, 1)  Measured angular velocity (rad/s)

noise = gyro.sigG/sqrt(dt)*randn(3, 1);

wGyro = (eye(3) + gyro.S)*wTrue + gyro.bias + noise;

end