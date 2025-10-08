import networkx as nx
import matplotlib.pyplot as plt
import argparse
import random
import itertools
import os
import pandas as pd
import numpy as np
from networkx.algorithms import community
from scipy import stats

# calculates neighborhood ovewrlap for plotting line thickness
def neighborhood_overlap(G):
    overlap = {}
    for u, v in G.edges():
        neighbors_u = set(G.neighbors(u)) - {v}
        neighbors_v = set(G.neighbors(v)) - {u}
        intersection = neighbors_u & neighbors_v
        union = neighbors_u | neighbors_v
        if len(union) == 0:
            overlap[(u, v)] = 0  # avoid division by zero
        else:
            overlap[(u, v)] = len(intersection) / len(union)
    return overlap

# returns the sign of the edge given
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

# reusable function for simualting failures and used in robustness check
def simulate_failures(G, k):
    Gcopy = G.copy()
    edges = list(Gcopy.edges())

    # removes a random set of k edges
    if k >= len(edges):
        to_remove = edges
    else:
        to_remove = random.sample(edges, k)  # uses random to decide which to remove
    Gcopy.remove_edges_from(to_remove)

    # finds the largest components before and after to determine average shortest path
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
    num_comp = nx.number_connected_components(Gcopy)

    # betweenness centrality impact by computing change per node, and averaging it
    bc_before = nx.betweenness_centrality(G)
    bc_after = nx.betweenness_centrality(Gcopy)
    delta_bc = {node: bc_after[node] - bc_before[node] for node in G.nodes()}
    avg_change = sum(abs(v) for v in delta_bc.values()) / len(delta_bc)

    # find how often the components persist
    # closer to 1 means more persisting and towards 0 means disrupted
    # Use greedy modularity communities
    original_communities = list(community.greedy_modularity_communities(G))
    new_communities = list(community.greedy_modularity_communities(Gcopy))

    # quick macro for calculations
    def jaccard_similarity(set1, set2):
        return len(set1 & set2) / len(set1 | set2)
    
    # calculate average similarity
    count = 0
    total = 0
    for i, orig in enumerate(original_communities):
        best_match = max(new_communities, key=lambda new: jaccard_similarity(orig, new))
        similarity = jaccard_similarity(orig, best_match)
        total += similarity
        count += 1
    
    cluster_persist = total / count # find average

    # find smallest and largest component
    components = list(nx.connected_components(Gcopy))

    # Largest component by number of nodes
    largest_cc = max(components, key=len)

    # Smallest component by number of nodes
    smallest_cc = min(components, key=len)

    return orig_aspl, new_aspl, num_comp, avg_change, len(smallest_cc), len(largest_cc), cluster_persist


def main(argv= None):

    # All the arg parser setup to read the input gml and command line inputs.
    # Some are optional while others exit the program if left out.
    p = argparse.ArgumentParser(description="Graph Analysis CLI (GML + NetworkX)")
    p.add_argument("graph_file", help="Input graph in .gml format")
    p.add_argument("--components", default=None, help="Partition graph into N components using Girvan-Newman")
    p.add_argument("--plot", choices=['C', 'N', 'P'], help="Plot mode: C=clustering, N=neighborhood overlap, P=attributes")
    p.add_argument("--verify_homophily", action='store_true', help="Perform homophily t-test using node attribute 'color'")
    p.add_argument("--verify_balanced_graph", action='store_true', help="Verify signed graph balance (edge attribute 'sign')")
    p.add_argument("--output", help="Save final graph to this .gml file")
    p.add_argument("--simulate_failures", default=None, help="Simulate removal of k random edges")
    p.add_argument("--robustness_check", nargs='?', const='1', help="Perform robustness check: provide k or leave blank to use 1")
    p.add_argument("--temporal_simulation", help="CSV file of edge changes for temporal simulation")
    p.add_argument("--split_arg_dir", help="Splitting the gml files by component")
    args = p.parse_args(argv)

    # load the graph from input .gml
    # graph will be an undirected graph object
    G = nx.read_gml(args.graph_file)

    # graph to modify if the components flag is called and you want to save it to an output
    output_graph = nx.Graph()

    # variables for creating temporal csv
    temporal = np.array() # stores our actions for later to convert into csv
    timestamp = 1 # iterable timestamp tracking each action

    # Ensure undirected for many algorithms unless the graph loaded is directed
    if isinstance(G, nx.DiGraph):
        print("Notice: Loaded directed graph, converting to undirected for analysis.")
        G = nx.Graph(G)

    # simulate_failures flag checked
    if args.simulate_failures:
        print('\n-----Simulate Failures-----')
        try:
            k = int(args.simulate_failures)
        except:
            print('Please give a valid integer as the k value for simulate_failures. Exiting program.')
            exit()
        results = simulate_failures(G, k)
        
        print(f'original avg shortest path: {results[0]}')
        print(f'new avg shortest_path: {results[1]}')
        print(f'change in avg shortest path: {None if (results[0] is None or results[1] is None) else (results[1] - results[0])}')
        print(f'number of disconnected components: {results[2]}')
        print(f'avg change in betweenness of nodes: {results[3]}')
    
    # robustness_check flag checked
    if args.robustness_check:
        print('\n-----Robustness Check-----')

        # simulates failures 100 times and reports average values using the given k value
        try:
            k = int(args.robustness_check)
        except:
            print('Please give a valid integer as the k value for robustness_check. Exiting program.')
            exit()
        
        # runs the simulation failure 100 times
        connected_components = 0
        min_component_size = None
        max_component_size = None
        persist = 0
        for i in range(100):
            results = simulate_failures(G, k)
            connected_components += results[2]

            if min_component_size == None:
                min_component_size = results[4]
            else:
                if min_component_size > results[4]:
                    min_component_size = results[4]

            if max_component_size == None:
                max_component_size = results[5]
            else:
                if max_component_size < results[5]:
                    max_component_size = results[5]

            persist += results[6]
        
        # average and print
        connected_components /= 100
        persist /= 100

        print(f'average number of components: {connected_components}')
        print(f'smallest component: {min_component_size}')
        print(f'largest component: {max_component_size}')
        print(f'component persistence(closer to 1 means persistent, to 0 is disrupted): {persist}')
        
        
    if args.components:
        n = args.components
        while True:
            try:
                n = int(n)
                if n<= 0:
                    print("Components n must be greater than 0, please input your new value")
                    n = input()
                else:
                    break
            except:
                print("please enter a number for Components n")
                n = input()
        
        comp_gen = nx.community.girvan_newman(G)

        # Iterate until we have n components
        for communities in itertools.islice(comp_gen, n - 1):
            result = [set(c) for c in communities]
        
        output = nx.Graph()
        for nodes in result:
            subgraph = G.subgraph(nodes).copy()
            output.add_nodes_from(subgraph.nodes(data=True))
            output.add_edges_from(subgraph.edges(data=True))
        output_graph = output

        # Return the last partition reached
        if args.split_arg_dir:
            os.makedirs(args.split_arg_dir, exist_ok=True)

            for i, nodes in enumerate(result, start=1):
                subgraph = G.subgraph(nodes).copy()
                filename = os.path.join(args.split_arg_dir, f"component_{i}.gml")
                nx.write_gml(subgraph, filename)
                print(f"Exported component {i} → {filename}")
        

    if args.plot:
        # input validation, skips if value is not either C, N, or P
        if args.plot not in ['C', 'N', 'P']:
            print("Plot flag was not given either one of C, N, or P. Please give a correct input, skipping plotting.")
        else:
            plot_mode = args.plot
            pos = nx.spring_layout(G, seed=42)
            plt.figure(figsize=(10, 8))
            if plot_mode == 'C':
                # node sizes: clustering
                clustering = nx.clustering(G)
                deg = dict(G.degree())
                sizes = [300 + 2000 * clustering.get(n, 0) for n in G.nodes()]
                colors = [deg.get(n, 0) for n in G.nodes()]
                nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color=colors, cmap='viridis')
                nx.draw_networkx_edges(G, pos, alpha=0.6)
                nx.draw_networkx_labels(G, pos, font_size=8)
                plt.title("Clustering Coefficient (node size) and Degree (color)")
            elif plot_mode == 'N':
                # edge thickness proportional to neighborhood overlap
                no = neighborhood_overlap(G)
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
            plt.show()
            plt.close()

    if args.verify_homophily:
        isHomophily = False
        
        same = 0
        different = 0

        for u, v in G.edges():
            if 'color' in G.nodes[u] and 'color' in G.nodes[v]:
                if G.nodes[u]['color'] == G.nodes[v]['color']:
                    same += 1
                else:
                    different += 1

        # if the graph given has no color assignments to differentiate, skips the homophily test
        if same + different == 0:
            print("No color value in nodes to test for homophily, skipping homophily test.")
        else:
            observed_ratio = same / (same + different)
            
            component_ids = [G.nodes[n]['color'] for n in G.nodes()]
            random_ratios = []
            for _ in range(1000):  # 1000 random permutations
                shuffled = random.sample(component_ids, len(component_ids))
                shuffled_components = dict(zip(G.nodes(), shuffled))
                same_r = 0
                diff_r = 0
                for u, v in G.edges():
                    if shuffled_components[u] == shuffled_components[v]:
                        same_r += 1
                    else:
                        diff_r += 1
                random_ratios.append(same_r / (same_r + diff_r))
            mean_random = sum(random_ratios) / len(random_ratios)
            
            t_stat, p_value = stats.ttest_1samp(random_ratios, observed_ratio)
            if observed_ratio > mean_random and p_value < 0.05:
                isHomophily = True
            print("homophily test")
            print(isHomophily)

    if args.verify_balanced_graph:
        balanced = True
        color_map = {} 

        for node in G.nodes():
            if node not in color_map:
                color_map[node] = 0
                q = [node]

                while q:
                    u = q.pop(0)
                    for v in G.neighbors(u):
                        sign = G[u][v].get('sign', 1)  # default +1 if missing
                        if v not in color_map:
                            if sign > 0:
                                color_map[v] = color_map[u]  # same set
                            else:
                                color_map[v] = 1 - color_map[u]  # opposite set
                            q.append(v)
                        else:
                            # Check consistency
                            if sign > 0 and color_map[v] != color_map[u]:
                                balanced = False
                                break
                            if sign < 0 and color_map[v] == color_map[u]:
                                balanced = False
                                break
                    if not balanced:
                        break
            if not balanced:
                break
        print("balance test")
        print(balanced)
    
    if args.output:
        # adds metadata for component id's
        components = list(nx.connected_components(output_graph))
        for component_id, comp in enumerate(components):
            nx.set_node_attributes(output_graph, {n: component_id for n in comp}, name="component_id")

        # all nodes get indicated true or false for isolation
        isolated_nodes = list(nx.isolates(output_graph))
        nx.set_node_attributes(output_graph, "False", name="isolated")
        nx.set_node_attributes(output_graph, {n: "True" for n in isolated_nodes}, name="isolated")
        nx.write_gml(output_graph, f"{args.output}")

    # # Temporal simulation output to csv
    # if args.temporal_simulation:
        

if __name__ == "__main__":
    main()