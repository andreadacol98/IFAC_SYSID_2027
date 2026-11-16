function y_lim = fcn_ylim(ev, const)
%% fcn_ylim
%
% DESCRIPTION
%   Common y-limits for the figures of a set of runs, as in Fig. 4: from
%   the lowest value reached by lambda_min(P_{t-1}) in any run up to the
%   highest upper bound pi_hi (t+1)^{-beta}, attained at t = 1.
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
