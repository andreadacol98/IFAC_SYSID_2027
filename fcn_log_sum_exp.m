function log_s = fcn_log_sum_exp(log_a, log_b)
%% fcn_log_sum_exp
%
% DESCRIPTION
%   Numerically stable evaluation of log(exp(log_a) + exp(log_b)),
%   elementwise, i.e. the log of a sum of two nonnegative quantities given
%   only through their logs. It never forms exp(log_a) or exp(log_b)
%   directly, so sums of quantities that would overflow/underflow in linear
%   scale (e.g. the effective noise variances phi_t^i ~ q_t^{-alpha_i} of a
%   deep tree, see fcn_phi) stay representable in log scale. A zero term is
%   passed as log(0) = -Inf.
%
% INPUT
%   log_a   array in [-Inf, Inf), log of the first term
%   log_b   array in [-Inf, Inf), log of the second term (same size as
%           log_a, or compatible with it)
%
% OUTPUT
%   log_s   array, log(exp(log_a) + exp(log_b)) (-Inf where both terms are 0)
%

m = max(log_a, log_b);
log_s = m + log1p(exp(-abs(log_a - log_b)));
log_s(m == -Inf) = -Inf;    % both terms are exp(-Inf) = 0 (and -Inf + Inf = NaN above)

end
