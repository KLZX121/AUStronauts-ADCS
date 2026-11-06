function PlotEclipsePatch(eclipseData, t)
%PlotEclipsePatch Plots an eclipse shaded region to an existing plot
%   PlotEclipsePatch(eclipseData, t)
%
%   Use after the desired plot.
%
%   Inputs
%   eclipseData     (n, 1)  List of n timesteps of boolean eclipse flags
%   t               (n, 1)  List of timesteps

    d = diff([0; eclipseData; 0]);
    % first eclipse sample of each run
    iStart = find(d == 1);
    % first non-eclipse sample (clamped to last)
    iEnd   = min(find(d == -1), numel(t));

    if isempty(iStart)
        return
    end

    t = t(:);
    xS = t(iStart)';
    xE = t(iEnd)';
    patchX = [xS; xE; xE; xS]; % 4 x k, one column per patch

    hold on
    yl = ylim;
    ylim(yl);
    patchY = repmat([yl(1); yl(1); yl(2); yl(2)], 1, numel(iStart));

    patch(patchX, patchY, [0.5 0.5 0.5], 'FaceAlpha', 0.2, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off', 'PickableParts', ...
        'none', 'HitTest', 'off');

    hold off

end