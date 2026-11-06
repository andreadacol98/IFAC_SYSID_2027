function net = fcn_tree(E, epsilon)
%% fcn_tree
%
% DESCRIPTION
%   Builds the directed tree (V, E) of Sec. 3.1 and its augmented graph
%   (V_bar, E_bar). The agents are V = {1, ..., r}; an edge (j, i) in E
%   means that j is a child of i, so every edge points toward the root r,
%   which has no parent, while every other agent has exactly one. The
%   fictitious agent 0, the signal s_t(x), is a child of every leaf:
%
%     E_bar = E u {(0, i): i in L}.
%
%   For each agent i, also computes its children C_bar^i in the augmented
%   graph (0 for a leaf) and its shortest distance to a leaf
%
%     delta_i = 0                              if i in L,
%     delta_i = 1 + min_{j in C_bar^i} delta_j otherwise,
%
%   which sets the exponent alpha_i of Lemma 1, and an ordering of V in
%   which every agent comes after all its children: a round of
%   word-of-mouth estimation proceeds along it, from the leaves to the root
%   (Fig. 2).
%
% INPUT
%   E         |E| x 2, edges (j, i) of the tree, labels in {1, ..., r}, root r
%             (empty for a single agent, r = 1)
%   epsilon   scalar, observation noise variance epsilon^{ij} > 0 of every
%             edge (j, i) in E_bar (Assumption 2)
%
% OUTPUT
%   net   struct with fields
%           r         number of agents, = label of the root
%           E         |E| x 2, edges (j, i)
%           L         leaves (column)
%           C         r x 1 cell, C{i} = C_bar^i (column), = 0 for a leaf
%           epsilon   r x 1 cell, epsilon{i}(k) = epsilon^{ij}, j = C{i}(k)
%           delta     r x 1, delta_i
%           order     r x 1, agents sorted so that children come before parents
%

if isempty(E)
    E = zeros(0, 2);
end

%% VALIDATION (Sec. 3.1: a directed tree rooted at r, every other agent with exactly one parent)
if size(E, 2) ~= 2 || any(E(:) < 1) || any(E(:) ~= round(E(:)))
    error("fcn_tree:invalidEdges", "E must be a |E| x 2 array of positive integer labels");
end
r = max([1; E(:)]);
if ~isequal(sort(E(:, 1)).', 1:r-1)
    error("fcn_tree:notTree", ...
        "every agent but the root r = %d must be the child in exactly one edge, and the root in none", r);
end
if any(E(:, 1) == E(:, 2))
    error("fcn_tree:selfLoop", "E contains a self-loop");
end
if ~(isscalar(epsilon) && epsilon > 0 && isfinite(epsilon))
    error("fcn_tree:invalidNoise", "epsilon must be a finite, strictly positive scalar");
end

%% CHILDREN AND LEAVES (augmented graph: the signal, agent 0, is the only child of each leaf)
C = cell(r, 1);
for i = 1:r
    C{i} = E(E(:, 2) == i, 1);
end
L = find(cellfun(@isempty, C));
for i = L.'
    C{i} = 0;
end

%% ORDER FROM THE LEAVES TO THE ROOT (reversed breadth-first visit from r)
order = r;
k = 1;
while k <= numel(order)
    if C{order(k)}(1) ~= 0
        order = [order; C{order(k)}];
    end
    k = k + 1;
end
if numel(order) ~= r
    error("fcn_tree:notConnected", "not every agent is connected to the root r = %d (E has a cycle)", r);
end
order = flipud(order);

%% SHORTEST DISTANCE TO A LEAF
delta = zeros(r, 1);
for i = order.'
    if C{i}(1) ~= 0
        delta(i) = 1 + min(delta(C{i}));
    end
end

%% NETWORK STRUCT
net.r = r;
net.E = E;
net.L = L;
net.C = C;
net.epsilon = cellfun(@(C_i) epsilon * ones(numel(C_i), 1), C, "UniformOutput", false);
net.delta = delta;
net.order = order;

end
