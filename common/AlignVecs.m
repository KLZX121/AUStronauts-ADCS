function R = AlignVecs(u, v)
%R = AlignVecs(u, v)
%
%   Returns the rotation matrix required to rotate one vector to another
%   Uses Rodrigues' rotation formula
%
%   Inputs
%   u   (3, 1)  Vector to align to
%   v   (3, 1)  Vector to rotate
%
%   Outputs
%   R   (3, 3)  Rotation matrix where u = Rv

    u = u./norm(u);
    v = v./norm(v);

    n = cross(v, u);
    n = n./norm(n);
    c = dot(v, u);

    if norm(n) < 1e-10
        if c > 0
            % parallel
            R = eye(3);
        else
            % anti parallel
            [~, i] = min(abs(v));
            e = zeros(3, 1);
            e(i) = 1;

            n = cross(v, e);
            n = n./norm(n);

            R = 2*(n*n') - eye(3);
        end
    else
        % general case
        R = eye(3) + Skew(n) + Skew(n)^2/(1+c);
    end
end