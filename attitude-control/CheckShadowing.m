function isShadowed = CheckShadowing(v, srcPanel, targetSurfaces)
%isShadowed = CheckShadowing(v, srcPanel, targetSurfaces)
%   Checks whether a panel is shadowed by other panels
%   
%   Inputs
%   v               (3, 1)      Outward direction of disturbance source (velocity or sun)
%   srcPanel        (struct)    Panel to check for shadowing
%   targetSurfaces  (struct[])  List of panels to check against for shadowing
%
%   Outputs
%   isShadowed      (bool)      Whether the panel is shadowed

isShadowed = false;
for j = 1:length(targetSurfaces)
    targetPanel = targetSurfaces(j);

    c = dot(v, targetPanel.normal);

    if abs(c) < eps
        % perpendicular (no intersection
        continue
    end
    P1 = targetPanel.vertices(1, :)';
    P2 = targetPanel.vertices(2, :)';
    P4 = targetPanel.vertices(4, :)';

    t = dot((P1 - srcPanel.centroid), targetPanel.normal) / c;

    if t <= 0
        % target is behind source
        continue
    end

    % intersection point
    r = srcPanel.centroid + v*t;

    % check if within bounds
    a = P2 - P1;
    b = P4 - P1;
    q = r - P1;

    u = dot(q, a) / dot(a, a);
    w = dot(q, b) / dot(b, b);

    isShadowed = u >= 0 && u <= 1 && w >= 0 && w <= 1;

    if isShadowed
        return
    end
end

end