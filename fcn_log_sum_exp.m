function log_s = fcn_log_sum_exp(log_a, log_b)
%% fcn_log_sum_exp
%
% DESCRIPTION
%   Numerically stable log(exp(log_a) + exp(log_b)), elementwise: the log
%   of a sum of two nonnegative quantities given through their logs. A zero
%   term is passed as log(0) = -Inf.
%
% INPUT
%   log_a   array in [-Inf, Inf), log of the first term
%   log_b   array in [-Inf, Inf), log of the second term (size compatible with log_a)
%
% OUTPUT
%   log_s   array, log(exp(log_a) + exp(log_b)) (-Inf where both terms are 0)
%

m = max(log_a, log_b);
log_s = m + log1p(exp(-abs(log_a - log_b)));
log_s(m == -Inf) = -Inf;    % both terms are 0

end
