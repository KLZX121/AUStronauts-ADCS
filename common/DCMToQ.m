function q = DCMToQ(A)
%DCMToQ Converts a direct-cosine-matrix (DCM) into a quaternion
%   q = DCMToQ(A)
%
%   The formula used is from Markley 2014 pg.48, rewritten for scalar-first
%   quaternions
%
%   Inputs
%   A   (3, 3)  Direct cosine matrix
%
%   Outputs
%   q   (4, 1)  Scalar-first quaternion

[~, i] = max([trace(A); A(1, 1); A(2, 2); A(3, 3)]);

if (i == 1)
    q = [
        1 + trace(A);
        A(2, 3) - A(3, 2);
        A(3, 1) - A(1, 3);
        A(1, 2) - A(2, 1);
    ];
elseif (i == 2)
    q = [
        A(2, 3) - A(3, 2);
        1 + 2*A(1, 1) - trace(A);
        A(1, 2) + A(2, 1);
        A(1, 3) + A(3, 1);
    ];
elseif (i == 3)
    q = [
        A(3, 1) - A(1, 3);
        A(2, 1) + A(1, 2);
        1 + 2*A(2, 2) - trace(A);
        A(2, 3) + A(3, 2);
    ];
elseif (i == 4)
    q = [
        A(1, 2) - A(2, 1);
        A(3, 1) + A(1, 3);
        A(3, 2) + A(2, 3);
        1 + 2*A(3, 3) - trace(A);
    ];
end

q = q./norm(q);

end