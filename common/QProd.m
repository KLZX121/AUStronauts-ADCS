function q = QProd(q1, q2)
%QProd Calculates the Shuster product of two quaternions
%   q = QProd(q1, q2)
%
%   This function produces the quaternion:
%   q = q1*q2
%   Let q2 be a rotation from A -> B, and q1 be a rotation from B -> C.
%   The resulting q is then a rotation from A -> C
%   Note that this product differs to the Hamilton product.
%   Src: Markley 2014 pg. 37 (rewritten for scalar-first)
%
%   Inputs
%   q1  (4, 1)  First quaternion
%   q2  (4, 1)  Second quaternion
%
%   Outputs
%   q   (4, 1)  Quaternion product

q = [
    q1(1)*q2(1) - dot(q1(2:4), q2(2:4));
    q2(1)*q1(2:4) + q1(1)*q2(2:4) - cross(q1(2:4), q2(2:4));
];

end