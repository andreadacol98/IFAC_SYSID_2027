function ev = fcn_eig(traj, const)
%% fcn_eig
%
% DESCRIPTION
%   Extracts from one run the smallest and the largest eigenvalue of the
%   common error covariance P_{t-1}, t = 1, ..., T+1 (from P_0 to P_T), to
%   be compared with the bounds of Theorem 1, and the range of
%   phi_t^i q_t^{alpha_i} over t = 1, ..., T for every agent, to be
%   compared with the constants [L_i, U_i] of Lemma 1 (fcn_constants).
%
% INPUT
%   traj    traj struct, after fcn_run; uses P, q, phi
%   const   struct, explicit constants (fcn_constants); uses alpha
%
% OUTPUT
%   ev      struct with fields
%             t                1 x (T+1), times t = 1, ..., T+1
%             lambda           2 x (T+1), [lambda_min; lambda_max] of P_{t-1}
%             phi_lo, phi_hi   r x 1, min and max of phi_t^i q_t^{alpha_i} over t = 1, ..., T
%

%% EIGENVALUES OF P_{t-1}
n_t = size(traj.P, 3);
ev.t = 1:n_t;
if size(traj.P, 1) == 2
    % closed form for symmetric 2 x 2 matrices, all t at once: (a + d)/2 -+ sqrt(((a - d)/2)^2 + b^2)
    a = reshape(traj.P(1, 1, :), 1, n_t);
    b = reshape(traj.P(1, 2, :), 1, n_t);
    d = reshape(traj.P(2, 2, :), 1, n_t);
    ev.lambda = (a + d)/2 + [-1; 1] .* hypot((a - d)/2, b);
else
    ev.lambda = NaN(2, n_t);
    for t = ev.t
        e = eig(traj.P(:, :, t));
        ev.lambda(:, t) = [min(e); max(e)];
    end
end

%% RANGE OF phi_t^i q_t^{alpha_i} (Lemma 1)
phi_scaled = traj.phi(:, 2:end) .* traj.q(2:end).^const.alpha;
ev.phi_lo = min(phi_scaled, [], 2);
ev.phi_hi = max(phi_scaled, [], 2);

end
