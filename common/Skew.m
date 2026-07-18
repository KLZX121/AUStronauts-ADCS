function M = Skew(v)
%Skew Generates the cross product matrix of a vector
%   M = Skew(v)
%
%   Inputs
%   v   (3, 1)  Vector
%
%   Outputs
%   M   (3, 3)  Skew symmetric matrix

M = [0      -v(3)   v(2);
     v(3)   0       -v(1);
     -v(2)  v(1)     0;
];

end