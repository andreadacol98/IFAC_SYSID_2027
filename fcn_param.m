function param = fcn_param()
%% fcn_param
%
% DESCRIPTION
%   Problem parameters of the numerical examples (Sec. 6). The unknown
%   x = [x_1; x_2], with mean x_hat0 and covariance P0 > 0, is revealed
%   through the scalar signal (1)
%
%     s_t(x) = c_t' x,   c_t = [cos(omega t); rho sin(omega t)],
%
%   with rho > 0 and omega not in pi Z (otherwise x_2 is never observed).
%
%   Assumption 1 holds in closed form for every t: with W = t_bar + 1,
%
%     c_lo = min(1, rho^2) (W - |sin(W omega) / sin(omega)|) / 2,
%     c_hi = max(1, rho^2) (W + |sin(W omega) / sin(omega)|) / 2,
%
%   where W is the smallest window with |sin(W omega) / sin(omega)| <= W/2.
%
% INPUT
%   none
%
% OUTPUT
%   param   struct with fields
%             n               dimension of x (= 2)
%             x_hat0          n x 1, mean of x
%             P0              n x n, covariance of x (P0 > 0)
%             T               simulation length
%             omega           frequency of the regressor (not in pi Z)
%             rho             gain of the sine channel of the regressor (> 0)
%             c               n x (T+1), regressor c_t (column t+1 <-> time t, NaN at t = 0)
%             t_bar           Assumption 1 holds over windows of t_bar + 1 steps
%             c_lo, c_hi      constants of Assumption 1
%             c_sup           bound on ||c_t||_2^2, = max(1, rho^2)
%             q_max           bound on c_t' P0 c_t >= q_t, = lambda_max(D P0 D), D = diag(1, rho)
%             epsilon         observation noise variance epsilon^{ij}, the same on every edge of E_bar
%             E_A             |E| x 2, edges (j, i) of tree A, delta_r = 2 (Fig. 3(a))
%             E_B             |E| x 2, edges (j, i) of tree B, delta_r = 1 (Fig. 3(b))
%             tree            3 x 1 string, tree of each run (panels (a)-(c) of Fig. 4)
%             eta             3 x 1, normalization parameter eta of (4) of each run
%

%% HORIZON
param.T = 1e5;

%% REGRESSOR (1), c_t = [cos(omega t); rho sin(omega t)]
param.n = 2;
param.omega = 1;
param.rho = 2;
t = 1:param.T;
param.c = [NaN(param.n, 1), [cos(param.omega * t); param.rho * sin(param.omega * t)]];

%% PRIOR (x in R^n, mean x_hat0, covariance P0)
param.x_hat0 = zeros(param.n, 1);
param.P0 = diag([2, 1]);

%% STRONG PERSISTENCY OF EXCITATION (Assumption 1)
s_omega = abs(sin(param.omega));
if s_omega < 1e-8
    error("fcn_param:notPersistentlyExciting", ...
        "omega = %g is a multiple of pi: sin(omega t) = 0 for every t, so x_2 is never observed (Assumption 1)", param.omega);
end
W = 2:min(ceil(2 / s_omega), param.T);                  % candidate window lengths t_bar + 1
dirichlet = abs(sin(W * param.omega)) / s_omega;        % |sum of exp(2 i omega tau)| over a window of W steps
W = W(find(dirichlet <= W / 2, 1));                     % smallest admissible window
if isempty(W)
    error("fcn_param:shortHorizon", ...
        "omega = %g is too close to a multiple of pi to excite every direction within T = %d steps", param.omega, param.T);
end
param.t_bar = W - 1;
dirichlet = abs(sin(W * param.omega)) / s_omega;
D = diag([1, param.rho]);                               % c_t = D [cos(omega t); sin(omega t)]
param.c_lo = min(1, param.rho^2) * (W - dirichlet) / 2;
param.c_hi = max(1, param.rho^2) * (W + dirichlet) / 2;
param.c_sup = max(1, param.rho^2);                      % ||c_t||_2^2 <= c_sup
param.q_max = max(eig(D * param.P0 * D));               % c_t' P0 c_t <= q_max

%% NOISE (Assumption 2)
param.epsilon = 1;       % variance epsilon^{ij} of e_t^{ij}, every edge (j, i) of E_bar

%% NETWORKS (edges (j, i): j is a child of i; the root r has the largest label, Sec. 3.1)
% tree A, delta_r = 2 (Fig. 3(a))
%
%          7
%        /   \
%       2     6
%       |     |
%       1     5
%            / \
%           3   4
%
param.E_A = [1 2; 2 7; 3 5; 4 5; 5 6; 6 7];

% tree B, delta_r = 1 (Fig. 3(b))
%
%          6
%        /   \
%       1     5
%             |
%             4
%            / \
%           2   3
%
param.E_B = [1 6; 2 4; 3 4; 4 5; 5 6];

%% RUNS (panels (a)-(c) of Fig. 4)
param.tree = ["A"; "A"; "B"];
param.eta = [0; 1; 1];   % beta = 1/11, 1/7, 1/3 (Table 1)

%% VALIDATION (Sec. 3: P0 > 0, c_t ~= 0, Assumption 1, Assumption 2: epsilon > 0, (4): eta in [0, 1])
if ~(param.rho > 0) || isinf(param.rho)
    error("fcn_param:invalidGain", "rho must be finite and strictly positive, got %g", param.rho);
end
[~, flag] = chol(param.P0);
if ~isequal(param.P0, param.P0.') || flag ~= 0
    error("fcn_param:invalidP0", "P0 must be symmetric positive definite");
end
if any(sum(param.c(:, 2:end).^2, 1) == 0)
    error("fcn_param:zeroRegressor", "c_t must be nonzero for every t");
end
if ~(param.c_lo > 0)
    error("fcn_param:notPersistentlyExciting", ...
        "c_t is not strongly persistently exciting over windows of %d steps (c_lo = %g)", param.t_bar + 1, param.c_lo);
end
[c_lo_T, c_hi_T] = fcn_pe(param.c(:, 2:end), param.t_bar);   % over the windows in the horizon
tol = 1e-9;
if c_lo_T < param.c_lo * (1 - tol) || c_hi_T > param.c_hi * (1 + tol) || ...
        max(sum(param.c(:, 2:end).^2, 1)) > param.c_sup * (1 + tol) || ...
        max(sum(param.c(:, 2:end) .* (param.P0 * param.c(:, 2:end)), 1)) > param.q_max * (1 + tol)
    error("fcn_param:closedFormMismatch", "the closed-form constants of Assumption 1 do not bound the regressor");
end
if ~(param.epsilon > 0) || isinf(param.epsilon)
    error("fcn_param:invalidNoise", "epsilon must be finite and strictly positive, got %g", param.epsilon);
end
if numel(param.tree) ~= numel(param.eta) || any(~(param.eta >= 0 & param.eta <= 1))
    error("fcn_param:invalidEta", "every run needs one tree and one eta in [0, 1]");
end

end

function [c_lo, c_hi] = fcn_pe(c, t_bar)
% Extreme eigenvalues of the sums of c_t c_t' over every window
% {t, ..., t + t_bar} that fits in c (n x T).
[n, T] = size(c);
if t_bar + 1 > T
    error("fcn_param:shortHorizon", "T = %d is shorter than a window of t_bar + 1 = %d steps", T, t_bar + 1);
end
K = T - t_bar;                                                  % number of windows
S = zeros(n, n, K);                                             % S(:, :, t) = sum_{tau = t}^{t + t_bar} c_tau c_tau'
for a = 1:n
    for b = a:n
        cs = [0, cumsum(c(a, :) .* c(b, :))];
        S(a, b, :) = reshape(cs(t_bar + 2:end) - cs(1:K), 1, 1, K);
        S(b, a, :) = S(a, b, :);
    end
end
if n == 2
    % closed form for symmetric 2 x 2 matrices
    a = reshape(S(1, 1, :), 1, K);
    b = reshape(S(1, 2, :), 1, K);
    d = reshape(S(2, 2, :), 1, K);
    lambda = (a + d)/2 + [-1; 1] .* hypot((a - d)/2, b);
else
    lambda = zeros(2, K);
    for k = 1:K
        ev = eig(S(:, :, k));
        lambda(:, k) = [min(ev); max(ev)];
    end
end
c_lo = min(lambda(1, :));
c_hi = max(lambda(2, :));
end
