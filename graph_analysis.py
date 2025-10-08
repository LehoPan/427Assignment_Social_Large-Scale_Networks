import networkx as nx
import matplotlib.pyplot as plt
import argparse
import itertools

def main(argv= None):

    # All the arg parser setup to read the input gml and command line inputs.
    # Some are optional while others exit the program if left out.
    p = argparse.ArgumentParser(description="Graph Analysis CLI (GML + NetworkX)")
    p.add_argument("graph_file", help="Input graph in .gml format")
    p.add_argument("--components", type=int, default=None, help="Partition graph into N components using Girvan-Newman")
    p.add_argument("--plot", choices=['C', 'N', 'P'], help="Plot mode: C=clustering, N=neighborhood overlap, P=attributes")
    p.add_argument("--verify_homophily", action='store_true', help="Perform homophily t-test using node attribute 'color'")
    p.add_argument("--verify_balanced_graph", action='store_true', help="Verify signed graph balance (edge attribute 'sign')")
    p.add_argument("--output", help="Save final graph to this .gml file")
    p.add_argument("--simulate_failures", type=int, default=None, help="Simulate removal of k random edges")
    p.add_argument("--robustness_check", nargs='?', const='1', help="Perform robustness check: provide k or leave blank to use 1")
    p.add_argument("--temporal_simulation", help="CSV file of edge changes for temporal simulation")
    args = p.parse_args(argv)

    # main graph will be an undirected graph object
    graph = nx.Graph()
        
    
    

if __name__ == "__main__":
    main()