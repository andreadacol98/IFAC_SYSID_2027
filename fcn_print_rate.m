function fcn_print_rate(name, net, eta, const, ev)
%% fcn_print_rate
%
% DESCRIPTION
%   Prints, for one run, the rate exponent of Theorem 1, the explicit
%   constants of Lemma 1 (root) and Theorem 1 (fcn_constants), and two
%   checks over the simulated horizon: every agent satisfies
%   L_i <= phi_t^i q_t^{alpha_i} <= U_i, and at every t
%   pi_lo (t+1)^{-beta} <= lambda_min(P_{t-1}) <= lambda_max(P_{t-1}) <= pi_hi (t+1)^{-beta}.
%
% INPUT
%   name    string, name of the tree
%   net     struct, directed tree (fcn_tree); uses r, delta
%   eta     normalization parameter of the run
%   const   struct, explicit constants (fcn_constants)
%   ev      struct, eigenvalues of the run (fcn_eig)
%

r = net.r;
tol = 1e-9;                                             % relative rounding tolerance of the checks
scale = (ev.t + 1).^const.beta;                         % (t+1)^beta
lemma_ok = all(ev.phi_lo >= const.L * (1 - tol) & ev.phi_hi <= const.U * (1 + tol));
theorem_ok = all(scale .* ev.lambda(1, :) >= const.pi_lo * (1 - tol)) && ...
    all(scale .* ev.lambda(2, :) <= const.pi_hi * (1 + tol));

fprintf("tree %s, delta_r = %d, eta = %g: alpha_r = %g, beta = %.4g\n", name, net.delta(r), eta, const.alpha(r), const.beta);
fprintf("  Theorem 1: [pi_lo, pi_hi] = [%.3g, %.3g]\n", const.pi_lo, const.pi_hi);
fprintf("  Lemma 1, root: [L_r, U_r] = [%.3g, %.3g], phi_t^r q_t^alpha_r in [%.3g, %.3g] over the run, limit %.3g\n", ...
    const.L(r), const.U(r), ev.phi_lo(r), ev.phi_hi(r), const.phi_inf(r));
fprintf("  checks over the run: Lemma 1 holds for every agent: %s, Theorem 1 holds at every t: %s\n", ...
    mat2str(lemma_ok), mat2str(theorem_ok));

end
