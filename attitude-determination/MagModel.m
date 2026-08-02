function bMag = MagModel(x, d)
%MagModel Magnetometer Model
%   bMag = MagModel(x, d)
%   
%   Inputs
%   x               (:, 1)      state vector
%   d               (struct)    data struct
%   .jD             (1, 1)      Julian Date of epoch
%   .iR             (1, 3)      indices of ECI pos in state vector
%   .iQ             (1, 4)      indices of attitude quaternion
%   .mag.sigma      (3, 1)      1-sigma noise for each axis
%   .mag.bias       (3, 1)      biases
%   .mag.D          (3, 3)      scale factor and non-orthogonality matrix
%   .mag.O          (3, 3)      DCM for sensor to body frame rotation
%   
%   Outputs
%   bMag  (3, 1)  vector magnetic field measurement (T) (body frame)

A = QToDCM(x(d.iQ));

bRef = IGRFECI(x(d.iR), d.jD);

noise = d.mag.sigma.*randn(3, 1);

bMag = (eye(3) + d.mag.D) \ (d.mag.O'*A*bRef + d.mag.bias + noise);

end
