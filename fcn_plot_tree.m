function fig = fcn_plot_tree(net, plot_param)
%% fcn_plot_tree
%
% DESCRIPTION
%   Draws the augmented tree (V_bar, E_bar) as in Fig. 1: the fictitious
%   agent 0, i.e. the signal s_t(x), at the bottom, the root r at the top,
%   and every edge pointing from a child to its parent. Agent 0 is drawn as
%   a square and the root as a diamond.
%
% INPUT
%   net          struct, directed tree (fcn_tree); uses r, E, L, delta
%   plot_param   struct, plotting style (fcn_plot_param)
%
% OUTPUT
%   fig          handle of the created figure
%

tail = [zeros(numel(net.L), 1); net.E(:, 1)];   % children, 0 = signal
head = [net.L; net.E(:, 2)];                    % parents
G = digraph(tail + 1, head + 1);                % node k + 1 <-> agent k

fig = figure("Visible", "on");                  % separate figure window, not embedded in the live script
ax = axes(fig);
h = plot(ax, G, "Layout", "layered", "Direction", "up", ...
    "NodeLabel", cellstr(string(0:net.r)), "NodeColor", "k", "EdgeColor", [0.4, 0.4, 0.4], ...
    "MarkerSize", 6, "LineWidth", 1, "ArrowSize", 8, "NodeFontSize", plot_param.font_size);
highlight(h, 1, "Marker", "s", "MarkerSize", 8);            % agent 0, the signal
highlight(h, net.r + 1, "Marker", "d", "MarkerSize", 8);    % root

axis(ax, "off");
title(ax, sprintf("$\\delta_r = %d$ (agent 0: signal $s_t(\\mathbf{x})$, root $r = %d$)", net.delta(net.r), net.r), ...
    "Interpreter", plot_param.interpreter, "FontSize", plot_param.font_size);

end
