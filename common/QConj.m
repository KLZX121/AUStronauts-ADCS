function qC = QConj(q)
%QConj Calculates the conjugate of a quaternion
%   qC = QConj(q)
%
%   The conjugate of a quaternion is equal to its inverse if it is a unit
%   quaternion
%
%   Inputs
%   q    (4, 1)  Quaternion
%
%   Outputs
%   qC   (4, 1)  Conjugate quaternion

qC = zeros(4, 1);

qC(1) = q(1);
qC(2:4) = -q(2:4);

end