function plot_param = fcn_plot_param(n_runs)
%% fcn_plot_param
%
% DESCRIPTION
%   Plotting style shared by all figures: one color and one marker per run
%   (cycled if n_runs exceeds the palette), line width, fonts and axes options.
%
% INPUT
%   n_runs       number of runs to be told apart
%
% OUTPUT
%   plot_param   struct with fields
%                  color         n_runs x 3, RGB color of each run ("gem" palette)
%                  marker        n_runs x 1 string, marker of each run
%                  marker_size   marker size
%                  line_width    line width
%                  font_size     font size of labels, legends and ticks
%                  interpreter   text interpreter ("latex")
%                  grid          axes grid ("on"/"off")
%

%% RUN COLORS
palette = orderedcolors("gem");
plot_param.color = palette(mod(0:n_runs-1, size(palette, 1)) + 1, :);

%% RUN MARKERS
plot_param.marker = repmat("o", n_runs, 1);
plot_param.marker_size = 6;

%% LINES
plot_param.line_width = 1.5;

%% FONTS
plot_param.font_size = 10;
plot_param.interpreter = "latex";

%% AXES
plot_param.grid = "on";

end
