function y_lim = fcn_ylim(ev, const)
%% fcn_ylim
%
% DESCRIPTION
%   Common y-limits for the figures of a set of runs, so that they have
%   the same size and can be compared: from the lowest value ever reached
%   by the smallest eigenvalue lambda_min(P_{t-1}) of any run, up to the
%   highest value of the upper bounds pi_hi (t+1)^{-beta}, attained at
%   t = 1.
%
% INPUT
%   ev      cell array, eigenvalues of the runs (fcn_eig)
%   const   cell array of the same size, explicit constants of the runs (fcn_constants)
%
% OUTPUT
%   y_lim   1 x 2, [lowest lambda_min, highest upper bound]
%

y_lim = [min(cellfun(@(e) min(e.lambda(1, :)), ev(:))), ...
    max(cellfun(@(c) c.pi_hi * 2^(-c.beta), const(:)))];

end
