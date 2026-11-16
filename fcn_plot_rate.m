function [fig, ax] = fcn_plot_rate(ev, const, plot_param, k, y_lim)
%% fcn_plot_rate
%
% DESCRIPTION
%   One panel of Fig. 4, in its own figure window: on log-log axes, the
%   smallest and the largest eigenvalue of P_{t-1}, t = 1, ..., T+1,
%   against the upper bound of Theorem 1,
%
%     lambda_min(P_{t-1}) <= lambda_max(P_{t-1}) <= pi_hi (t+1)^{-beta},
%
%   with the region below the bound shaded. The lower bound
%   pi_lo (t+1)^{-beta} lies many decades below the eigenvalues and is not
%   drawn (it is printed by fcn_print_rate). Markers are placed at 12
%   log-spaced times, and the bound is drawn through 200 log-spaced times.
%
% INPUT
%   ev           struct, eigenvalues of the run (fcn_eig)
%   const        struct, explicit constants of the run (fcn_constants)
%   plot_param   struct, plotting style (fcn_plot_param)
%   k            index of the color and marker of the run in plot_param
%   y_lim        (optional) 1 x 2, y-limits, e.g. common to several runs
%                (fcn_ylim); default: from the lowest lambda_min of the run
%                to the highest value of its upper bound
%
% OUTPUT
%   fig          handle of the created figure
%   ax           handle of its axes
%

t = ev.t;
y_lo = ev.lambda(1, :);                                 % lambda_min(P_{t-1})
y_hi = ev.lambda(2, :);                                 % lambda_max(P_{t-1})
t_b = unique([round(logspace(0, log10(t(end)), 200)), t(end)]);   % times at which the bound is drawn
b_hi = const.pi_hi * (t_b + 1).^(-const.beta);         % upper bound of Theorem 1
color = plot_param.color(k, :);
marker = plot_param.marker(k);
idx = unique(round(logspace(0, log10(t(end)), 12)));    % markers at log-spaced times
if nargin < 5
    y_lim = [min(y_lo), max(b_hi)];
end

fig = figure("Visible", "on");                          % separate figure window
ax = axes(fig);
hold(ax, "on");
ax.XScale = "log";
ax.YScale = "log";

% region below the upper bound, down to the bottom of the axes (drawn past it and clipped)
fill(ax, [t_b, fliplr(t_b)], [b_hi, y_lim(1) / 10 * ones(size(t_b))], color, "FaceAlpha", 0.15, "EdgeColor", "none");
h_b = plot(ax, t_b, b_hi, "--", "Color", color, "LineWidth", plot_param.line_width);

% eigenvalues at every t: largest with filled markers, smallest with empty markers
h_hi = plot(ax, t, y_hi, "Color", color, "LineWidth", plot_param.line_width, ...
    "Marker", marker, "MarkerIndices", idx, "MarkerSize", plot_param.marker_size, "MarkerFaceColor", color);
h_lo = plot(ax, t, y_lo, "Color", color, "LineWidth", plot_param.line_width, ...
    "Marker", marker, "MarkerIndices", idx, "MarkerSize", plot_param.marker_size, "MarkerFaceColor", "w");

legend(ax, [h_hi, h_lo, h_b], ...
    ["$\lambda_{\max}(\mathbf{P}_{t-1})$", "$\lambda_{\min}(\mathbf{P}_{t-1})$", "$\bar{\pi}\,(t+1)^{-\beta}$"], ...
    "Interpreter", plot_param.interpreter, "FontSize", plot_param.font_size, "Location", "northeast");

xlim(ax, [t(1), t(end)]);
ylim(ax, y_lim);
xlabel(ax, "$t$", "Interpreter", plot_param.interpreter, "FontSize", plot_param.font_size);
box(ax, "on");
grid(ax, plot_param.grid);
ax.TickLabelInterpreter = plot_param.interpreter;
hold(ax, "off");

end
