function traj = fcn_traj(param, net)
%% fcn_traj
%
% DESCRIPTION
%   Initializes one trajectory struct: draws the unknown x and the noise
%   sequences, and allocates the quantities filled by fcn_run. The noise
%   e_t^{ij} corrupting the observations of agent i is drawn upfront, one
%   row per child j in C_bar^i. Column k of every time-indexed field
%   corresponds to time k-1, t = 0, ..., T. x and the noises are Gaussian,
%   while the draft only needs their first two moments (Assumption 2).
%
% INPUT
%   param   struct, problem parameters (fcn_param); uses n, x_hat0, P0, T
%   net     struct, directed tree (fcn_tree); uses r, C, epsilon
%
% OUTPUT
%   traj    struct with fields
%             x       n x 1, unknown x ~ N(x_hat0, P0)
%             e       r x 1 cell, e{i}(k, t+1) = e_t^{ij}, j = net.C{i}(k),
%                     ~ N(0, epsilon^{ij}) (NaN at t = 0)
%             x_hat   n x (T+1), common estimate x_hat_t broadcast by the root (column 1: x_hat0)
%             P       n x n x (T+1), common error covariance P_t (page 1: P0)
%             q       1 x (T+1), q_t = c_t' P_{t-1} c_t (NaN at t = 0)
%             phi     r x (T+1), effective noise variance phi_t^i of every agent (NaN at t = 0)
%

%% UNKNOWN (x ~ N(x_hat0, P0))
traj.x = param.x_hat0 + chol(param.P0, "lower") * randn(param.n, 1);

%% NOISE SEQUENCES (exogenous, generated upfront, one row per edge (j, i) in E_bar)
traj.e = cell(net.r, 1);
for i = 1:net.r
    n_C = numel(net.C{i});
    traj.e{i} = [NaN(n_C, 1), sqrt(net.epsilon{i}) .* randn(n_C, param.T)];
end

%% COMMON ESTIMATE (broadcast by the root, shared by all agents)
traj.x_hat = NaN(param.n, param.T+1);
traj.P = NaN(param.n, param.n, param.T+1);
traj.x_hat(:, 1) = param.x_hat0;
traj.P(:, :, 1) = param.P0;

%% DETERMINISTIC QUANTITIES (Proposition 1 (ii), filled by fcn_run)
traj.q = NaN(1, param.T+1);
traj.phi = NaN(net.r, param.T+1);

end
