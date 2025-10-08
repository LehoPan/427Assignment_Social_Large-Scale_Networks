# 427Assignment Social Large-Scale Networks
By Leo Pan 030025552 and Nicolas Piker 029966545

## Dependencies
We installed and utilized the python `networkx`, `scipy`, and `matplotlib` packages, in addition to the base python packages `argparse`, `random`, `os`, and `itertools`.

## Implementations
For this project we implemented following functions/arguments:
- `--components`
    - Input validation checks if a validn integer is given, otherwise skips
    - We utilized the `nx.community.girvan_newman()` function to remove the most used edges until desired number of components is reached.
    - Side note that if the robustness_check or simulate_failures flag is set they will both run before components.
- `--split_output_dir`
    - Used enumerate to create subgraphs and write each component to their own .gml file.
- `--plot`
    - Validates input is one of "C", "P", or "N". If not it will skip the plotting.
    - Uses `matplotlib` to plot each node and edge in the graph. Computes the clustering with `nx.clustering()`, and neighborhood overlap with a custom function instead of `networkx`.
- `--verify_homophily`
    - Only works if the nodes in the given .gml file have a `color` attribute, and groups them by that. If there is no color attribute found, test is skipped, otherwise calculates homphily for the graph. 
- `--verify_balanced_graph`
    - Checks for a `sign` value for each edge in the given .gml. If none is found defaults to a strong `1` value, `-1` means weak bond. 
- `--output`
    - Outputs the graph and writes component ids to each node, and if they are isolated. Important if the `--component` flag is flipped, it'll correctly output the new graph that has been divided into more components. Exports to the file with .gml appended.
- `--simulate_failures`
    - Since it is used by robustness check, it has been made its own function for ease of reusability. Randomly removes k number of edges with the `random` library. and calculates shortest paths, number of components, and average changing of betweenness.
    - Also returns smallest and largest components, and how well the clusters persist, as robustness check needs that additional information.
- `--robustness_check`
    - Calls simulated failures function 100 times with the given k value.
    - Keeps track of the smallest and largest components across all the iterations.
    - Averages how well the components persist when determined with jaccard similarity. It will be between 0 and 1, and the closer it is to 1 means more persisting components, and 0 means more disrupted components by the edge failures.
- `--temporal_simulation`
    - Creates a csv file with the edges from the input graph and timestamps each edge in the order that they are added in the gml file.

## Running the program
You can utilize any number of the flags with an input graph.gml file. This graph will need a `color` attribute for every node if you want to check for homophily, the naming can be anything you want. A `sign` attribute on every edge, of either `1` or `-1` if you want to check for balanced graphs. We noticed that there is no `--plot T` requirement, but the temporal flag is still there. So we did not implement a T argument for plot, but we still have an output .csv generated with the temporal_simulation flag.

e.g. commands:
`python ./graph_analysis.py graph.gml --components 3 --plot C --simulate_failures 5 --output output.gml`

`python ./graph_analysis.py graph.gml --verify_homophily --verify_balanced_graph --output output.gml`