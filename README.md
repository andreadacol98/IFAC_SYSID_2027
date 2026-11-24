# Do Hierarchical Networks Actually Learn? Insights from Word-of-Mouth Estimation

## Abstract
We study recursive estimation of a constant state vector by a hierarchical network of asymmetrically informed linear estimators, modeled as a directed tree.
Only the leaves observe the regression signal, whereas every other estimator learns by *word of mouth* from noisy signal estimates relayed by its children.
Authority within the hierarchy is expressed through common prior information: after each round, every estimator adopts the estimate and error covariance broadcast by the root.
The more accurate this common information, the larger the effective noise variance on local observations, which makes estimation intrinsically self-defeating.
The resulting performance is the subject of this paper.
We show that recursive least-squares (RLS) estimators do learn the state in mean square, albeit slowly.
Under strong persistency of excitation, the broadcast error covariance decays polynomially in time, at a rate that we characterize by upper and lower bounds of matching order.
The exponent is determined by the smallest number of relays between a leaf and the root, and it decreases exponentially with this number.

## Additional Information
This repository provides commented MATLAB code for simulating Word-of-Mouth estimation in hierarchical networks, reproducing the numerical examples of the paper, and checking the theoretical bounds on the common error covariance $\mathbf{P}_t$ over a finite simulation horizon.

For each configuration, the rate exponent $\beta$ and the explicit constants of Lemma 1 and Theorem 1 are computed before the simulation. The code then checks the corresponding inequalities at every applicable simulated time step, up to a relative tolerance for numerical rounding. These finite-horizon checks assess the stated bounds; they do not, by themselves, independently estimate an asymptotic convergence exponent.

The effective observation-noise variances can scale as $q_t^{-\alpha_i}$, where the exponents $\alpha_i$ increase exponentially with the shortest distance $\delta_i$ from an agent to a leaf. Selected intermediate calculations are performed in the log domain to mitigate overflow and underflow; this does not eliminate all floating-point limitations.

The implementation also generates a realization of the unknown state and observation-noise sequences to simulate the estimation protocol. Under the model, however, the common error covariance is deterministic, so its theoretical trajectory and the covariance plots do not depend on the random seed. The seed fixes the simulated stochastic trajectory.

The code accompanies the paper submitted to the **21st IFAC Symposium on System Identification (SYSID 2027)**:
[🔗 Do Hierarchical Networks Actually Learn? Insights from Word-of-Mouth Estimation][paper].

## Requirements
**MATLAB R2024b** or later (tested on R2024b). No specific toolboxes are needed.

## Installation
Download the repository as a ZIP archive, or clone it using:

```bash
git clone https://github.com/andreadacol98/IFAC_SYSID_2027.git
```

## File Overview
Run the live script `main.mlx` to reproduce the numerical examples. No interactive input is required. The default configurations are defined in `fcn_param.m`, which can be edited to change the regressor parameters, observation-noise variance, simulation horizon, normalization parameters, or network topologies.

The repository contains the following **MATLAB** files:

| File | Description |
|---|---|
| `main.mlx` | Main live script for the numerical examples and their plots. |
| `fcn_agent.mlx` | Local update of agent $i$ (Fig. 2): inverse-variance weighting (4), RLS update (6a)–(6b), and transmission of the signal estimate to the parent. |
| `fcn_constants.mlx` | Computes the exponents $\alpha_i$ and $\beta$, the agent-wise constants of Lemma 1, and the constants of Theorem 1. |
| `fcn_eig.m` | Extracts the extreme eigenvalues of $\mathbf{P}_{t-1}$ and the simulated ranges of $\varphi_t^i q_t^{\alpha_i}$ for each agent. The latter products are evaluated directly in floating-point arithmetic. |
| `fcn_log_sum_exp.m` | Evaluates $\log(e^a+e^b)$ in a numerically stable way for log-domain calculations. |
| `fcn_param.m` | Defines the model and simulation parameters, computes explicit excitation bounds for the sinusoidal regressor, and checks the main parameter constraints and model assumptions. |
| `fcn_phi.mlx` | Computes the parameters of Proposition 1 ($a_t^{ij}$, $b_t^{ij}$, $\varphi_t^{ij}$, $\varphi_t^i$, $a_t^i$) from the leaves to the root, using log-domain calculations where needed. |
| `fcn_plot_param.m` | Sets the shared plotting style: colors for different runs, circular markers, line widths, fonts, and axes options. |
| `fcn_plot_rate.m` | Plots the smallest and largest eigenvalues of $\mathbf{P}_{t-1}$ together with the upper bound of Theorem 1 on log-log axes, in a separate figure window. |
| `fcn_plot_tree.m` | Draws the augmented directed tree with the fictitious signal agent 0 at the bottom, the root at the top, and edges directed toward the root. |
| `fcn_print_rate.m` | Prints the exponent $\beta$, the constants of Theorem 1, and the root's Lemma 1 constants; checks the Lemma 1 bounds for every agent and the Theorem 1 bounds throughout the simulated horizon. |
| `fcn_run.mlx` | Runs the Word-of-Mouth estimation protocol up to time $T$, updating agents from the leaves to the root and broadcasting the common pair $(\hat{\mathbf{x}}_t,\mathbf{P}_t)$ after each round. |
| `fcn_traj.m` | Samples the unknown state and observation-noise sequences, initializes the common estimate and covariance, and allocates the trajectory arrays. |
| `fcn_tree.m` | Builds and validates a directed tree from its edges, treats the fictitious signal agent as a child of each leaf, and determines children, leaves, shortest distances $\delta_i$, and an update order from children to parents. |
| `fcn_ylim.m` | Computes common vertical plotting limits so that covariance plots can be compared on the same scale. |

## Usage
1. Add the repository folder to your **MATLAB** path.
2. Open and run `main.mlx`.

The default parameters in `fcn_param.m` correspond to Section 6 of the paper: $\mathbf{x}\in\mathbb{R}^2$, prior mean $\hat{\mathbf{x}}_0=\mathbf{0}$ and covariance $\mathbf{P}_0=\mathrm{diag}(2,1)$, regressor

```math
\mathbf{c}_t = [\cos(\omega t),\;\rho\sin(\omega t)]^\top,
\qquad \omega=1,\quad \rho=2,
```

unit observation-noise variance on each edge, and horizon $T=10^5$. The random seed is fixed by `rng(0)` for reproducibility of the sampled state and noise sequences; the common covariance trajectory is independent of that seed.

Two directed trees are used. Their roots are shown at the top, and edges are understood to point upward from children to parents. Agent 0 (the signal) and its edges to the leaves are omitted from these schematic drawings:

```text
      delta_r = 2                 delta_r = 1
       Fig. 3(a)                   Fig. 3(b)

           7                          6
         /   \                      /   \
        2     6                    1     5
        |     |                          |
        1     5                          4
             / \                        / \
            3   4                      2   3
```

The MATLAB network plots produced by `fcn_plot_tree.m` instead show the *augmented* trees, including the fictitious agent 0, in an upward layered layout. They therefore differ in presentation from the unaugmented network diagrams in Fig. 3 of the paper.

The experiments use the following settings:

- **Experiment 1 (normalization):** the network with $\delta_r=2$ is run with $\eta=0$ and $\eta=1$, giving $\beta=1/11$ and $\beta=1/7$, respectively.
- **Experiment 2 (depth):** the networks with $\delta_r=2$ and $\delta_r=1$ are compared at $\eta=1$, giving $\beta=1/7$ and $\beta=1/3$, respectively. The $\delta_r=2$, $\eta=1$ run is reused from Experiment 1.

The outputs correspond to the following paper figures. They are identified here by content rather than by MATLAB window number, which depends on the order in which `main.mlx` creates the figures:

| Output | Configuration | Paper |
|---|---|---|
| Network diagram | $\delta_r=2$ | Fig. 3(a) |
| Network diagram | $\delta_r=1$ | Fig. 3(b) |
| Covariance plot | $\delta_r=2$, $\eta=0$, $\beta=1/11$ | Fig. 4(a) |
| Covariance plot | $\delta_r=2$, $\eta=1$, $\beta=1/7$ | Fig. 4(b) |
| Covariance plot with the common y-axis limits of Experiment 2 | $\delta_r=2$, $\eta=1$, $\beta=1/7$ (same run) | Fig. 4(b) |
| Covariance plot | $\delta_r=1$, $\eta=1$, $\beta=1/3$ | Fig. 4(c) |

Each plotting function opens a separate MATLAB figure window. The covariance plots show the smallest and largest eigenvalues of $\mathbf{P}_{t-1}$ alongside the upper bound from Theorem 1. The lower bound is checked and its constant is printed but is not plotted, since it lies many orders of magnitude below the eigenvalues for the default settings. The plotted trajectories need not exhibit their asymptotic slopes over the finite horizon because convergence can be slow.

Increasing `param.T` in `fcn_param.m` from $10^5$ to $10^6$ extends the time axis by one decade and increases computational cost substantially (approximately proportionally to the horizon for the default settings).

[paper]: https://arxiv.org/abs/2611.20806
