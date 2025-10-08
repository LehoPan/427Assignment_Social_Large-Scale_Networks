import argparse
import os
import random
import math
from collections import deque, defaultdict
from typing import List, Tuple

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

try:
    import pandas as pd
except Exception:
    pd = None

try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False

# ---------------------------
# Utility functions
# ---------------------------

def load_graph(path: str) -> nx.Graph:
    """Load a GML graph (NetworkX supports node/edge attributes)."""
    print(f"Loading graph from '{path}'...")
    G = nx.read_gml(path)
    # Ensure undirected for many algorithms unless the graph loaded is directed
    if isinstance(G, nx.DiGraph):
        print("Warning: Loaded directed graph; converting to undirected for analysis.")
        G = nx.Graph(G)
    print(f"Graph loaded: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")
    return G

def write_graph(G: nx.Graph, path: str):
    print(f"Writing graph to '{path}'...")
    nx.write_gml(G, path)
    print("Write complete.")

# ---------------------------
# Metrics
# ---------------------------

def compute_clustering(G: nx.Graph, attr_name: str = "clustering") -> None:
    """Compute node clustering coefficients and attach as node attributes."""
    print("Computing clustering coefficients...")
    clustering = nx.clustering(G)
    nx.set_node_attributes(G, clustering, attr_name)
    avg = sum(clustering.values()) / len(clustering) if clustering else 0.0
    print(f"Avg clustering coefficient: {avg:.4f}")

def compute_neighborhood_overlap(G: nx.Graph, attr_name: str = "neighborhood_overlap") -> None:
    """
    Neighborhood overlap defined as:
      |N(u) ∩ N(v)| / |N(u) ∪ N(v)|
    for each edge (u,v). Stored as edge attribute.
    """
    print("Computing neighborhood overlap for edges...")
    for u, v in G.edges():
        Nu = set(G.neighbors(u)) - {v}
        Nv = set(G.neighbors(v)) - {u}
        inter = Nu & Nv
        union = Nu | Nv
        overlap = (len(inter) / len(union)) if union else 0.0
        G[u][v][attr_name] = overlap
    print("Neighborhood overlap computed.")

# ---------------------------
# Community (Girvan-Newman)
# ---------------------------

def girvan_newman_n_communities(G: nx.Graph, n: int) -> List[set]:
    """
    Partition graph into n communities using iterative Girvan-Newman:
    repeatedly remove edges of highest betweenness until number of
    connected components >= n, then return the components.
    """
    if n <= 0:
        raise ValueError("n must be > 0")
    if n == 1:
        return [set(G.nodes())]

    working = G.copy()
    components = list(nx.connected_components(working))
    iteration = 0
    while len(components) < n and working.number_of_edges() > 0:
        iteration += 1
        # compute edge betweenness on current graph
        eb = nx.edge_betweenness_centrality(working)
        # find max eb value and remove all edges with that value (ties)
        max_eb = max(eb.values())
        edges_to_remove = [e for e, val in eb.items() if math.isclose(val, max_eb)]
        print(f"[GN iter {iteration}] Removing {len(edges_to_remove)} edge(s) with betweenness {max_eb:.6f}")
        working.remove_edges_from(edges_to_remove)
        components = list(nx.connected_components(working))
    print(f"Girvan-Newman produced {len(components)} components.")
    return components

def label_communities(G: nx.Graph, communities: List[set], attr_name: str = "community"):
    mapping = {}
    for idx, comp in enumerate(communities):
        for node in comp:
            mapping[node] = idx
    nx.set_node_attributes(G, mapping, attr_name)
    print(f"Labeled {len(communities)} communities as node attribute '{attr_name}'.")

# ---------------------------
# Homophily test
# ---------------------------

def verify_homophily_ttest(G: nx.Graph, attribute: str = "color", sample_pairs: int = 10000) -> dict:
    """
    Verify homophily using a t-test:
       - For each existing edge, compute indicator I_edge = 1 if node attributes equal else 0
       - For random node pairs sample (not necessarily edges), compute I_random
       - Use t-test (independent samples) to check whether I_edge mean > I_random mean
    Returns dictionary with statistics.
    """
    if attribute not in next(iter(G.nodes(data=True)))[1]:
        raise KeyError(f"Attribute '{attribute}' not found on nodes. Node data sample: {list(G.nodes(data=True))[:3]}")

    print(f"Running homophily t-test on attribute '{attribute}'...")
    # list of attribute values per node
    nodes = list(G.nodes())
    attr_values = {n: G.nodes[n].get(attribute) for n in nodes}

    # edges indicators
    edge_inds = []
    for u, v in G.edges():
        edge_inds.append(1.0 if attr_values[u] == attr_values[v] else 0.0)

    # random pairs indicators
    rng = random.Random(42)
    random_inds = []
    max_pairs = sample_pairs
    n_nodes = len(nodes)
    for _ in range(min(max_pairs, sample_pairs)):
        a, b = rng.randrange(n_nodes), rng.randrange(n_nodes)
        if a == b:
            continue
        na, nb = nodes[a], nodes[b]
        random_inds.append(1.0 if attr_values[na] == attr_values[nb] else 0.0)
    # ensure same sizes for t-test convenience
    m = min(len(edge_inds), len(random_inds))
    edge_sample = edge_inds[:m]
    rand_sample = random_inds[:m]

    res = {}
    res['edge_mean'] = float(np.mean(edge_sample)) if edge_sample else 0.0
    res['random_mean'] = float(np.mean(rand_sample)) if rand_sample else 0.0
    # perform t-test if scipy available, otherwise do simple z-like approximate test
    if SCIPY_AVAILABLE:
        tstat, pval = stats.ttest_ind(edge_sample, rand_sample, equal_var=False)
        res['t_stat'] = float(tstat)
        res['p_value'] = float(pval)
        print(f"T-test result: t={tstat:.4f}, p={pval:.6f}")
    else:
        # fallback: compute difference and approximate z-score
        se = math.sqrt(np.var(edge_sample, ddof=1)/len(edge_sample) + np.var(rand_sample, ddof=1)/len(rand_sample))
        diff = res['edge_mean'] - res['random_mean']
        z = diff / se if se > 0 else float('nan')
        res['t_stat'] = z
        res['p_value'] = None
        print("scipy not available: returning approximate z-statistic (p_value=None).")
    return res

# ---------------------------
# Verify signed graph balance
# ---------------------------

def _edge_sign_value(G: nx.Graph, u, v):
    """Return sign as +1 or -1. Accept strings like '+' / '-' or numeric."""
    s = G[u][v].get('sign', None)
    if s is None:
        # default to positive if not specified
        return 1
    if isinstance(s, (int, float)):
        return 1 if s >= 0 else -1
    if isinstance(s, str):
        return 1 if s.strip() in ('1', '+', 'pos', 'positive') else -1
    return 1

def is_signed_graph_balanced(G: nx.Graph) -> Tuple[bool, dict]:
    """
    BFS-based check for structural balance:
    Attempt to label nodes group = +1 or -1 such that for every edge (u,v) with sign s:
      s = group[u] * group[v]  (i.e., s>0 => group[u]==group[v], s<0 => group[u]!=group[v])
    Returns (balanced_bool, groups_dict)
    """
    groups = {}
    for comp_nodes in nx.connected_components(G):
        for node in comp_nodes:
            if node in groups:
                continue
            # BFS from node
            groups[node] = 1
            q = deque([node])
            while q:
                u = q.popleft()
                for v in G.neighbors(u):
                    s = _edge_sign_value(G, u, v)
                    expected = groups[u] if s > 0 else -groups[u]
                    if v not in groups:
                        groups[v] = expected
                        q.append(v)
                    else:
                        if groups[v] != expected:
                            print(f"Conflict found on edge ({u},{v}) with sign {s}: assigned {groups[v]} vs expected {expected}.")
                            return False, groups
    return True, groups

# ---------------------------
# Simulate failures and robustness
# ---------------------------

def simulate_failures(G: nx.Graph, k: int, seed: int = None) -> dict:
    """
    Remove k random edges (without replacement) and compute:
      - change in average shortest path (increase/decrease)
      - number of disconnected components
      - impact on top betweenness centralities (list of top-k nodes betweenness)
    Returns dictionary with results and mutated graph copy under key 'G'.
    """
    if k <= 0:
        return {'G': G.copy(), 'info': 'k<=0, no removal'}

    rng = random.Random(seed)
    Gcopy = G.copy()
    edges = list(Gcopy.edges())
    if k >= len(edges):
        to_remove = edges
    else:
        to_remove = rng.sample(edges, k)
    Gcopy.remove_edges_from(to_remove)
    print(f"Removed {len(to_remove)} edges (k={k}).")

    # Average shortest path: only for largest connected component
    info = {}
    def avg_shortest_path_val(H):
        if H.number_of_nodes() <= 1:
            return 0.0
        if nx.is_connected(H):
            return nx.average_shortest_path_length(H)
        # compute for largest component
        largest = max(nx.connected_components(H), key=len)
        return nx.average_shortest_path_length(H.subgraph(largest))
    orig_aspl = avg_shortest_path_val(G)
    new_aspl = avg_shortest_path_val(Gcopy)
    info['orig_avg_shortest_path'] = orig_aspl
    info['new_avg_shortest_path'] = new_aspl
    info['delta_avg_shortest_path'] = None if (orig_aspl is None or new_aspl is None) else (new_aspl - orig_aspl)
    info['num_components'] = nx.number_connected_components(Gcopy)
    # betweenness centrality impact
    orig_bc = nx.betweenness_centrality(G)
    new_bc = nx.betweenness_centrality(Gcopy)
    # measure top changes
    def top_k_nodes(bc_dict, k=5):
        return sorted(bc_dict.items(), key=lambda kv: kv[1], reverse=True)[:k]
    info['orig_top_bc'] = top_k_nodes(orig_bc, 10)
    info['new_top_bc'] = top_k_nodes(new_bc, 10)
    info['removed_edges'] = to_remove
    return {'G': Gcopy, 'info': info}

def robustness_check(G: nx.Graph, k: int, trials: int = 20, seed: int = None) -> dict:
    """
    Perform multiple simulations of k random edge removals and report:
      - average number of connected components
      - max/min component sizes across trials
      - cluster persistence (how many nodes remain in the same community)
    """
    rng = random.Random(seed)
    components_stats = []
    component_sizes_all = []
    cluster_persistence_counts = []

    # baseline communities (single run)
    baseline_comps = [set(c) for c in nx.connected_components(G)]
    baseline_partition = {}
    for i, comp in enumerate(baseline_comps):
        for n in comp:
            baseline_partition[n] = i

    for t in range(trials):
        sim_seed = rng.randint(0, 2**30)
        res = simulate_failures(G, k, seed=sim_seed)
        Gsim = res['G']
        comps = [set(c) for c in nx.connected_components(Gsim)]
        components_stats.append(len(comps))
        sizes = [len(c) for c in comps]
        component_sizes_all.extend(sizes)
        # cluster persistence: how many nodes still in same baseline component id (heuristic)
        sim_partition = {}
        for i, comp in enumerate(comps):
            for n in comp:
                sim_partition[n] = i
        same = sum(1 for n in G.nodes() if baseline_partition.get(n) == sim_partition.get(n))
        cluster_persistence_counts.append(same / G.number_of_nodes())

    out = {
        'avg_num_components': float(np.mean(components_stats)),
        'max_component_size': int(max(component_sizes_all)) if component_sizes_all else 0,
        'min_component_size': int(min(component_sizes_all)) if component_sizes_all else 0,
        'avg_cluster_persistence': float(np.mean(cluster_persistence_counts))
    }
    print(f"Robustness over {trials} trials (k={k}): avg components {out['avg_num_components']:.2f}, "
          f"avg cluster persistence {out['avg_cluster_persistence']:.3f}")
    return out

# ---------------------------
# Plotting
# ---------------------------

def plot_graph(G: nx.Graph, plot_mode: str = "C", plot_attrs: dict = None, save_path: str = None):
    """
    plot_mode: 'C' clustering: node size = clustering, color = degree
               'N' neighborhood overlap: edge thickness = NO, color = sum degrees at endpoints
               'P' attributes: color nodes by attribute 'color' and edges by sign if present
    """
    pos = nx.spring_layout(G, seed=42)
    plt.figure(figsize=(10, 8))
    if plot_mode == 'C':
        # node sizes: clustering
        clustering = nx.get_node_attributes(G, 'clustering')
        deg = dict(G.degree())
        sizes = [300 + 2000 * clustering.get(n, 0) for n in G.nodes()]
        colors = [deg.get(n, 0) for n in G.nodes()]
        nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color=colors, cmap='viridis')
        nx.draw_networkx_edges(G, pos, alpha=0.6)
        nx.draw_networkx_labels(G, pos, font_size=8)
        plt.title("Clustering Coefficient (node size) and Degree (color)")
    elif plot_mode == 'N':
        # edge thickness proportional to neighborhood overlap
        no = nx.get_edge_attributes(G, 'neighborhood_overlap')
        # edge colors by sum of degrees
        edge_widths = [1 + 6 * no.get((u, v), no.get((v, u), 0.0)) for u, v in G.edges()]
        edge_colors = [G.degree(u) + G.degree(v) for u, v in G.edges()]
        nx.draw_networkx_nodes(G, pos, node_size=200)
        nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color=edge_colors, edge_cmap=plt.cm.plasma)
        nx.draw_networkx_labels(G, pos, font_size=8)
        plt.title("Neighborhood Overlap (edge thickness) and degree-sum (edge color)")
    elif plot_mode == 'P':
        # node color attribute 'color' (categorical) and edge sign show as dashed for negative
        node_colors_attr = nx.get_node_attributes(G, 'color')
        unique_vals = list(sorted(set(node_colors_attr.values()))) if node_colors_attr else []
        color_map = {}
        for i, val in enumerate(unique_vals):
            color_map[val] = i
        node_colors = [color_map.get(node_colors_attr.get(n, None), 0) for n in G.nodes()]
        # edges: positive solid, negative dashed
        pos_edges = [(u, v) for u, v in G.edges() if _edge_sign_value(G, u, v) > 0]
        neg_edges = [(u, v) for u, v in G.edges() if _edge_sign_value(G, u, v) < 0]
        nx.draw_networkx_nodes(G, pos, node_size=300, node_color=node_colors, cmap='tab10')
        nx.draw_networkx_edges(G, pos, edgelist=pos_edges)
        nx.draw_networkx_edges(G, pos, edgelist=neg_edges, style='dashed')
        nx.draw_networkx_labels(G, pos, font_size=8)
        plt.title("Attributes plot: node color by 'color', edges dashed if negative sign")
    else:
        print(f"Unknown plot mode '{plot_mode}' - no plot generated.")
        return
    plt.axis('off')
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    else:
        plt.show()
    plt.close()

# ---------------------------
# Temporal simulation & animation
# ---------------------------

def temporal_simulation_animate(G: nx.Graph, csv_file: str, output_gif: str = None, interval_ms: int = 500):
    """
    CSV format: source,target,timestamp,action
      - action: add or remove (case-insensitive)
    We'll animate changes over time and optionally save a GIF (requires imagemagick if saving via matplotlib).
    """
    if pd is None:
        raise RuntimeError("pandas required for temporal simulation. Please install pandas.")
    df = pd.read_csv(csv_file)
    # sort by timestamp
    df = df.sort_values(by='timestamp')
    events = df.to_records(index=False)
    # copy base graph
    Gt = G.copy()
    import matplotlib.animation as animation
    pos = nx.spring_layout(Gt, seed=42)  # static layout for consistency

    fig, ax = plt.subplots(figsize=(8, 6))

    def update_frame(i):
        ax.clear()
        if i >= len(events):
            return
        src, tgt, ts, action = events[i]
        action = str(action).strip().lower()
        if action in ("add", "insert", "create"):
            if not Gt.has_edge(src, tgt):
                Gt.add_edge(src, tgt)
        elif action in ("remove", "delete", "rm"):
            if Gt.has_edge(src, tgt):
                Gt.remove_edge(src, tgt)
        # draw current snapshot
        nx.draw_networkx_nodes(Gt, pos, ax=ax, node_size=150)
        nx.draw_networkx_edges(Gt, pos, ax=ax, alpha=0.8)
        ax.set_title(f"t={ts} action={action} ({src}-{tgt})")
        ax.axis('off')

    ani = animation.FuncAnimation(fig, update_frame, frames=len(events), interval=interval_ms, repeat=False)
    if output_gif:
        # saving requires imagemagick or pillow writer
        try:
            ani.save(output_gif, writer='pillow')
            print(f"Saved animation to {output_gif}")
        except Exception as e:
            print("Failed to save animation:", e)
    else:
        plt.show()

# ---------------------------
# Command-line interface
# ---------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Graph Analysis CLI (GML + NetworkX)")
    p.add_argument("graph_file", help="Input graph in .gml format")
    p.add_argument("--components", type=int, default=None, help="Partition graph into N components using Girvan-Newman")
    p.add_argument("--robustness_check", nargs='?', const='1', help="Perform robustness check: provide k or leave blank to use 1")
    p.add_argument("--simulate_failures", type=int, default=None, help="Simulate removal of k random edges")
    p.add_argument("--plot", choices=['C', 'N', 'P'], help="Plot mode: C=clustering, N=neighborhood overlap, P=attributes")
    p.add_argument("--verify_homophily", action='store_true', help="Perform homophily t-test using node attribute 'color'")
    p.add_argument("--verify_balanced_graph", action='store_true', help="Verify signed graph balance (edge attribute 'sign')")
    p.add_argument("--output", help="Save final graph to this .gml file")
    p.add_argument("--split_output_dir", help="If provided, export each component into separate .gml files in this dir")
    p.add_argument("--temporal_simulation", help="CSV file of edge changes for temporal simulation")
    p.add_argument("--plot_save", help="If set, save the plot to given path instead of showing")
    p.add_argument("--trials", type=int, default=20, help="Number of trials for robustness check")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("--components_attr_name", default="community", help="Node attribute name to store community id")
    return p.parse_args()

def main():
    args = parse_args()
    G = load_graph(args.graph_file)

    random.seed(args.seed)
    np.random.seed(args.seed)

    # compute baseline metrics
    compute_clustering(G)
    compute_neighborhood_overlap(G)

    # optionally simulate single failure removal immediately if --simulate_failures
    if args.simulate_failures:
        sim_res = simulate_failures(G, args.simulate_failures, seed=args.seed)
        G_after_fail = sim_res['G']
        print("Simulation info:", sim_res['info'])
        # Update G reference to the modified graph for subsequent actions
        G = G_after_fail

    # robust check: may be provided as --robustness_check k
    if args.robustness_check is not None:
        try:
            k = int(args.robustness_check)
        except Exception:
            k = 1
        rc = robustness_check(G, k, trials=args.trials, seed=args.seed)
        print("Robustness check summary:", rc)

    # Partition into N communities via Girvan-Newman
    if args.components is not None:
        n = args.components
        comps = girvan_newman_n_communities(G, n)
        label_communities(G, comps, attr_name=args.components_attr_name)
        # Optionally split output per component
        if args.split_output_dir:
            outdir = args.split_output_dir
            os.makedirs(outdir, exist_ok=True)
            for idx, comp in enumerate(comps):
                subG = G.subgraph(comp).copy()
                outpath = os.path.join(outdir, f"component_{idx}.gml")
                nx.write_gml(subG, outpath)
                print(f"Wrote component {idx} with {subG.number_of_nodes()} nodes to {outpath}")

    # Verify homophily
    if args.verify_homophily:
        try:
            res = verify_homophily_ttest(G, attribute='color')
            print("Homophily test result:", res)
        except Exception as e:
            print("Homophily test failed:", e)

    # Verify balance for signed graph
    if args.verify_balanced_graph:
        balanced, groups = is_signed_graph_balanced(G)
        print("Signed graph balanced?" , balanced)
        if balanced:
            nx.set_node_attributes(G, groups, 'signed_group')
            print("Assigned 'signed_group' node attributes for partitioning.")

    # Plot if requested
    if args.plot:
        plot_graph(G, plot_mode=args.plot, save_path=args.plot_save)

    # Temporal simulation (animation)
    if args.temporal_simulation:
        if pd is None:
            print("Temporal simulation requires pandas. Please pip install pandas.")
        else:
            # default gif path
            gif_out = os.path.splitext(args.temporal_simulation)[0] + "_anim.gif"
            try:
                temporal_simulation_animate(G, args.temporal_simulation, output_gif=gif_out)
            except Exception as e:
                print("Temporal animation failed:", e)

    # Final output save
    if args.output:
        # ensure updated node and edge attributes are present (clustering, neighborhood_overlap, community)
        write_graph(G, args.output)

    print("Analysis complete.")

if __name__ == "__main__":
    main()
