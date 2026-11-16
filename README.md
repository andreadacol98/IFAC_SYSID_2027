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
This repository provides commented code for simulating Word-of-Mouth estimation in hierarchical networks, and for verifying numerically the convergence rate of the common error covariance $\mathbf{P}_t$ broadcast by the root.
For every run, the rate exponent $\beta$ and the explicit constants of Lemma 1 and Theorem 1 are computed *before* the simulation, and the code then checks that the corresponding bounds hold at every simulated time step.
Since the effective noise variances grow as $q_t^{-\alpha_i}$, with $\alpha_i$ growing exponentially in the depth of the tree, these quantities are computed in log scale so that they remain within double precision.

The code was used to produce the numerical examples (Figures 3 and 4) of the following paper, submitted to the **IFAC Symposium on System Identification (SYSID 2027)**:
[🔗 Do Hierarchical Networks Actually Learn? Insights from Word-of-Mouth Estimation][paper].

## Requirements
**MATLAB R2024b** or later (tested on R2024b). No specific toolboxes are needed.

## Installation
Clone this repository by downloading the zip file. Alternatively, one can use:
```bash
git clone https://github.com/andreadacol98/IFAC_SYSID_2027.git
```

## File Overview
The code is organized so that the user only needs to run the live script `main.mlx`.
No input from the keyboard is required: all problem parameters are set in `fcn_param.m`, which can be edited to explore other settings (e.g., the frequency of the regressor, the noise variance, the horizon, or the trees).

Here is a high-level description of the **MATLAB** files in this repository:

| File | Description |
|---|---|
| `main.mlx` | Main live script to reproduce the numerical examples in the paper. |
| `fcn_agent.mlx` | Local update of agent $i$ (Fig. 2): inverse-variance weighting (4), RLS update (6a)–(6b), and the signal estimate transmitted to the parent. |
| `fcn_constants.mlx` | Computes, before any simulation, the exponents $\alpha_i$ and $\beta$, and the explicit constants of Lemma 1 (one pair per agent) and Theorem 1. |
| `fcn_eig.m` | Extracts the smallest and the largest eigenvalue of $\mathbf{P}_{t-1}$ from one run, and the range of $\varphi_t^i q_t^{\alpha_i}$ for every agent. |
| `fcn_log_sum_exp.m` | Numerically stable evaluation of $\log(e^a + e^b)$, used to keep quantities that would overflow or underflow representable in log scale. |
| `fcn_param.m` | Sets the problem parameters (regressor, prior, horizon, noise variance, normalization parameters $\eta$, and the two trees), computes the constants of Assumption 1 in closed form, and validates all inputs. |
| `fcn_phi.mlx` | Computes the parameters of Proposition 1 ($a_t^{ij}$, $b_t^{ij}$, $\varphi_t^{ij}$, $\varphi_t^i$, $a_t^i$) from the leaves to the root, in log scale. |
| `fcn_plot_param.m` | Sets the plotting style shared by all figures, with one color and one marker per run. |
| `fcn_plot_rate.m` | Plots the eigenvalues of $\mathbf{P}_{t-1}$ against the upper bound of Theorem 1 on log-log axes, for one run. |
| `fcn_plot_tree.m` | Draws the augmented tree, with the fictitious agent 0 (the signal) at the bottom and the root at the top. |
| `fcn_print_rate.m` | Prints the rate exponent and the explicit constants of one run, and checks that Lemma 1 and Theorem 1 hold at every simulated time step. |
| `fcn_run.mlx` | Runs Word-of-Mouth estimation up to time $T$: at every time step, the agents update from the leaves to the root, and the root broadcasts the common pair $(\hat{\mathbf{x}}_t, \mathbf{P}_t)$. |
| `fcn_traj.m` | Initializes one trajectory: draws the unknown state and the noise sequences, and allocates the quantities filled by `fcn_run.mlx`. |
| `fcn_tree.m` | Builds and validates a directed tree from its list of edges, and computes the children, the leaves, the distances $\delta_i$, and the processing order from the leaves to the root. |
| `fcn_ylim.m` | Computes common y-limits for the figures of a set of runs, so that they can be compared directly. |

## Usage
The same numerical examples of [🔗 Do Hierarchical Networks Actually Learn? Insights from Word-of-Mouth Estimation][paper] are considered, in order to show how the code is successfully run:
- Add the repository to your **MATLAB** path
- Open and run `main.mlx`
  - The default parameters in `fcn_param.m` are those of Section 6 of the paper: $\mathbf{x} \in \mathbb{R}^2$ with $\mathbf{P}_0 = \mathrm{diag}(2, 1)$, regressor $\mathbf{c}_t = [\cos(\omega t), \rho \sin(\omega t)]^\top$ with $\omega = 1$ and $\rho = 2$, unit noise variance on every edge, and horizon $T = 10^5$. The random seed is fixed by `rng(0)`.
  - Two directed trees are considered. Agent 0 (the signal) and its edges to the leaves are omitted below:

  ```text
     Tree A (delta_r = 1)       Tree B (delta_r = 2)
        Fig. 3(b)                  Fig. 3(a)

            6                          7
          /   \                      /   \
         1     5                    2     6
               |                    |     |
               4                    1     5
              / \                        / \
             2   3                      3   4
  ```

  - **Experiment 1 (normalization)** runs tree B with $\eta = 0$ and $\eta = 1$, for which $\beta = 1/11$ and $\beta = 1/7$, respectively.
  - **Experiment 2 (depth)** runs trees B and A with $\eta = 1$, for which $\beta = 1/7$ and $\beta = 1/3$, respectively. The run of tree B is reused from Experiment 1.
  - Every figure opens in its own window, not embedded in the live script, so that it can be edited with the figure tools. The lower bound of Theorem 1 is printed but not drawn, as it lies many orders of magnitude below the smallest eigenvalue.

  | Figure | Content | Paper |
  |---|---|---|
  | 1 | Tree A | Fig. 3(b) |
  | 2 | Tree B | Fig. 3(a) |
  | 3 | Tree B, $\eta = 0$, $\beta = 1/11$ | Fig. 4(a) |
  | 4 | Tree B, $\eta = 1$, $\beta = 1/7$ | Fig. 4(b) |
  | 5 | Tree B, $\eta = 1$, $\beta = 1/7$ (same run as Figure 4, on the axes of Experiment 2) | Fig. 4(b) |
  | 6 | Tree A, $\eta = 1$, $\beta = 1/3$ | Fig. 4(c) |

- The plots should be aligned with those provided in the link above.
- The rates of Theorem 1 are slow, so the horizon matters: increasing `param.T` in `fcn_param.m` to $10^6$ adds one decade to the plots, at roughly ten times the run time.

[paper]: https://arxiv.org/abs/XXXX.XXXXX
