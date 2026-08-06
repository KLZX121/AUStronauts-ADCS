classdef GyroModel < handle
properties
    bias (3, 1) double
    s (3, 1) double
    kU (3, 1) double
    kL (3, 1) double
    sigG (3, 1) double
    sigB (3, 1) double

    S (3, 3) double
end

methods(Static)
    %%% compose the calibration matrix S
    function S = SMatrix(s, kU, kL)
        S = [
            s(1) kU(1) kU(2);
            kL(1) s(2) kU(3);
            kL(2) kL(3) s(3)
            ];
    end
end
methods
    function o = GyroModel(calGyro)
        %GyroModel  Models a 3-axis rate gyro
        %   o = GyroModel(calGyro)
        %
        %   Inputs
        %   calGyro (struct)    Gyro Calibration parameters
        %   .bias   (3, 1)      Gyro bias (rad/s)
        %   .s      (3, 1)      Scale factors
        %   .kU     (3, 1)      Upper misalignments
        %   .kL     (3, 1)      Lower misalignments
        %   .S      (3, 3)      Scale factor/misalignment matrix
        %   .sigG   (3, 1)      1-sigma measurement noise (rad/s^(1/2))
        %   .sigB   (3, 1)      1-sigma bias walk (rad/s^(3/2))
        %
        %   Outputs
        %   o      (object)    GyroModel object

        o.bias = calGyro.bias;
        o.s = calGyro.s;
        o.kU = calGyro.kU;
        o.kL = calGyro.kL;
        o.sigG = calGyro.sigG;
        o.sigB = calGyro.sigB;

        o.S = o.SMatrix(calGyro.s, calGyro.kU, calGyro.kL);
    end

    function wGyro = Measurement(o, wTrue, dt)
        %Measurement Simulates a gyro measurement
        %   wGyro = Measurement(o, wTrue)
        %
        %   Inputs
        %   wTrue   (3, 1)  True angular velocity (rad/s)
        %   dt      (1, 1)  Timestep (s)
        %   
        %   Outputs
        %   wGyro   (3, 1)  Measured angular velocity (rad/s)

        noise = o.sigG./sqrt(dt).*randn(3, 1);

        wGyro = (eye(3) + o.S)*wTrue + o.bias + noise;
    end

    function o = PropagateBias(o, dt)
        %PropagateBias Propagates a random walk gyro bias
        %   o = PropagateBias(o, dt)
        %
        %   Inputs
        %   dt  (1, 1)      Timestep (s)
        %
        %   Outputs
        %   o   (object)    GyroModel object
    
        o.bias = o.bias + o.sigB.*sqrt(dt).*randn(3, 1);
    end
end
end

