function [fig, ax] = fcn_plot_rate(ev, const, net, eta, plot_param, k, y_lim)
%% fcn_plot_rate
%
% DESCRIPTION
%   Verifies the convergence rate of Theorem 1 for one run, on log-log
%   axes: the smallest and the largest eigenvalue of the common error
%   covariance P_{t-1} against the explicit upper bound of fcn_constants,
%
%     lambda_min(P_{t-1}) <= lambda_max(P_{t-1}) <= pi_hi (t+1)^{-beta},
%
%   with the region below the bound filled down to the bottom of the axes.
%   The figure opens in its own figure window, also when called from a
%   live script, so that it can be edited with the figure tools.
%   The bound is a straight line of slope -beta for t >> 1. The explicit
%   lower bound pi_lo (t+1)^{-beta} is not drawn: it lies many decades
%   below the eigenvalues (it is printed by fcn_print_rate). Every
%   eigenvalue of the sequence, t = 1, ..., T+1, is drawn; markers are
%   placed only at 12 times evenly spaced on the logarithmic time axis. The
%   bound is a smooth power law of t + 1, so it and the filled region are
%   drawn through 200 log-spaced times.
%
% INPUT
%   ev           struct, eigenvalues of the run (fcn_eig)
%   const        struct, explicit constants of the run (fcn_constants)
%   net          struct, directed tree of the run (fcn_tree); uses r, delta
%   eta          normalization parameter of the run
%   plot_param   struct, plotting style (fcn_plot_param)
%   k            index of the run in plot_param (color and marker)
%   y_lim        (optional) 1 x 2, y-limits, e.g. common to several runs
%                (fcn_ylim); default: from the lowest lambda_min of the run
%                to the highest value of its upper bound
%
% OUTPUT
%   fig          handle of the created figure
%   ax           handle of its axes (e.g. to share limits across figures)
%

t = ev.t;
y_lo = ev.lambda(1, :);                                 % lambda_min(P_{t-1})
y_hi = ev.lambda(2, :);                                 % lambda_max(P_{t-1})
t_b = unique([round(logspace(0, log10(t(end)), 200)), t(end)]);   % times at which the bound is drawn
b_hi = const.pi_hi * (t_b + 1).^(-const.beta);         % upper bound of Theorem 1
color = plot_param.color(k, :);
marker = plot_param.marker(k);
idx = unique(round(logspace(0, log10(t(end)), 12)));    % markers at log-spaced times
if nargin < 7
    y_lim = [min(y_lo), max(b_hi)];
end

fig = figure("Visible", "on");                          % separate figure window, not embedded in the live script
ax = axes(fig);
hold(ax, "on");
ax.XScale = "log";
ax.YScale = "log";

% region below the upper bound, down to the bottom of the axes (drawn past it and clipped)
fill(ax, [t_b, fliplr(t_b)], [b_hi, y_lim(1) / 10 * ones(size(t_b))], color, "FaceAlpha", 0.15, "EdgeColor", "none");
h_b = plot(ax, t_b, b_hi, "--", "Color", color, "LineWidth", plot_param.line_width);

% eigenvalues at every t, largest with filled markers, smallest with empty markers
h_hi = plot(ax, t, y_hi, "Color", color, "LineWidth", plot_param.line_width, ...
    "Marker", marker, "MarkerIndices", idx, "MarkerSize", plot_param.marker_size, "MarkerFaceColor", color);
h_lo = plot(ax, t, y_lo, "Color", color, "LineWidth", plot_param.line_width, ...
    "Marker", marker, "MarkerIndices", idx, "MarkerSize", plot_param.marker_size, "MarkerFaceColor", "w");

legend(ax, [h_hi, h_lo, h_b], ...
    ["$\lambda_{\max}(\mathbf{P}_{t-1})$", "$\lambda_{\min}(\mathbf{P}_{t-1})$", "$\bar{\pi}\,(t+1)^{-\beta}$"], ...
    "Interpreter", plot_param.interpreter, "FontSize", plot_param.font_size, "Location", "northeast");
title(ax, sprintf("$\\delta_r = %d,\\ \\eta = %g,\\ \\beta = %s$", net.delta(net.r), eta, fcn_rate_str(const.beta)), ...
    "Interpreter", plot_param.interpreter, "FontSize", plot_param.font_size);

xlim(ax, [t(1), t(end)]);
ylim(ax, y_lim);
xlabel(ax, "$t$", "Interpreter", plot_param.interpreter, "FontSize", plot_param.font_size);
ylabel(ax, "$\lambda(\mathbf{P}_{t-1})$", "Interpreter", plot_param.interpreter, "FontSize", plot_param.font_size);
box(ax, "on");
grid(ax, plot_param.grid);
ax.TickLabelInterpreter = plot_param.interpreter;
hold(ax, "off");

end

function s = fcn_rate_str(beta)
% beta as 1/k when 1/beta is an integer (e.g. eta in {0, 1}), with 3 significant digits otherwise.
if abs(1/beta - round(1/beta)) < 1e-9
    s = sprintf("1/%d", round(1/beta));
else
    s = sprintf("%.3g", beta);
end
end
