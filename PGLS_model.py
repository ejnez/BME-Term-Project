import pandas as pd
import numpy as np
from skbio import TreeNode
import statsmodels.api as sm


def brownian_covariance_matrix(tree, species_list):
    n = len(species_list)
    cov = np.zeros((n, n))
    root = tree

    for i, sp1 in enumerate(species_list):
        for j, sp2 in enumerate(species_list):
            node1 = tree.find(sp1)
            node2 = tree.find(sp2)
            mrca = tree.lowest_common_ancestor([node1, node2])
            cov[i, j] = mrca.distance(root)
    return cov


def run_pgls(order):
    print(f"\nProcessing order: {order}")
    
    # Load RES scores
    res_df = pd.read_csv(f"{order}_res_scores.csv", index_col=0)
    res_df.index = res_df.index.str.replace("_", " ")  # match tree labels
    
    # Extract lifespan from header
    lifespan = res_df.index.to_series().apply(lambda x: float(x.split("|")[1]))
    
    # Load tree
    tree = TreeNode.read(f"{order}_DBD_tree.nwk")
    
    # Reorder CSV to match tree tip order
    tree_tips = [tip.name for tip in tree.tips()]
    res_df = res_df.loc[tree_tips]
    lifespan = res_df.index.to_series().apply(lambda x: float(x.split("|")[1]))
    species_names = res_df.index.tolist()
    
    # Compute covariance matrix
    cov_matrix = brownian_covariance_matrix(tree, species_names)
    cov_matrix += np.eye(len(species_names)) * 1e-6  # numerical stability
    
   
    results = []
    for resi_col in res_df.columns:
        X = sm.add_constant(res_df[resi_col].values)
        y = lifespan.values
        gls_model = sm.GLS(y, X, sigma=cov_matrix)
        res = gls_model.fit()
        results.append({
            "Residue": resi_col,
            "Coefficient": res.params[1],
            "P_Value": res.pvalues[1],
            "T_Statistic": res.tvalues[1]
        })
    
    
    results_df = pd.DataFrame(results)
    output_file = f"{order}_PGLS_results.csv"
    results_df.to_csv(output_file, index=False)
    print(f"Saved results to {output_file}")


def main():
    orders = [
        "Artiodactyla",
        "Carnivora",
        "Chiroptera",
        "Primates",
        "Rodentia",
    ]
    
    for order in orders:
        run_pgls(order)


if __name__ == "__main__":
    main()
