from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt

def make_tree(input_file, output_file_prefix):
    """
    Build a phylogenetic tree from a trimmed DBD alignment.
    Uses Neighbor Joining (NJ) with identity distance.
    Saves both Newick and PNG.
    """

    # Load the alignment
    alignment = AlignIO.read(input_file, "fasta")

    # Identity distance
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Build NJ tree
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(distance_matrix)

    # Ladderize
    nj_tree.ladderize()

    # Save as Newick
    newick_file = f"{output_file_prefix}.nwk"
    Phylo.write(nj_tree, newick_file, "newick")
    print(f"Tree saved to {newick_file}")

    # --- Visualization ---
    fig = plt.figure(figsize=(20, 20))  # adjust size per preference
    ax = fig.add_subplot(1, 1, 1)

    # Draw tree
    Phylo.draw(nj_tree, axes=ax, do_show=False)

    # Reduce font size for tip labels
    for label in ax.get_ymajorticklabels():
        label.set_fontsize(9)

    png_file = f"{output_file_prefix}.png"
    plt.tight_layout()
    plt.savefig(png_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Tree figure saved to {png_file}")

    return nj_tree

def main():
    orders = [
        "Artiodactyla",
        "Carnivora",
        "Chiroptera",
        "Primates",
        "Rodentia",
    ]

    for order in orders:
        fasta_file = f"{order}_DBD_trimmed.fasta"
        output_prefix = f"{order}_DBD_tree"
        make_tree(fasta_file, output_prefix)

if __name__ == "__main__":
    main()
