function A = QToDCM(q)
%QToDCM Converts a quaternion into a direct-cosine-matrix (DCM)
%   A = QToDCM(q)
%
%   The formula used is from Markley 2014 pg.45, rewritten for scalar-first
%   quaternions in the form: q = [qs; qx; qy; qz]
%
%   Inputs
%   q   (4, 1)  Scalar-first quaternion
%
%   Outputs
%   A (3, 3)  Direct cosine matrix

qs = q(1);
qx = q(2);
qy = q(3);
qz = q(4);

A = [
  qx^2-qy^2-qz^2+qs^2   2*(qx*qy+qz*qs)         2*(qx*qz-qy*qs);
  2*(qy*qx-qz*qs)       -qx^2+qy^2-qz^2+qs^2    2*(qy*qz+qx*qs);
  2*(qz*qx+qy*qs)       2*(qz*qy-qx*qs)         -qx^2-qy^2+qz^2+qs^2;
];

end