import argparse
import random
from ete3 import Tree


def prune_random_leaves(treefile: str, outfile: str, max_leaves: int, seed: int = 42) -> None:
    """Prune random leaves from the tree until at most ``max_leaves`` remain."""
    tree = Tree(treefile, format=1)
    leaves = tree.get_leaves()
    if len(leaves) <= max_leaves:
        tree.write(outfile=outfile, format=1)
        return

    random.seed(seed)
    keep = random.sample(leaves, max_leaves)
    keep_names = [leaf.name for leaf in keep]
    tree.prune(keep_names, preserve_branch_length=True)
    tree.write(outfile=outfile, format=1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Randomly prune tree to limit number of leaves.")
    parser.add_argument("treefile", help="Input Newick tree file")
    parser.add_argument("outfile", help="Output Newick tree file")
    parser.add_argument("--max-leaves", type=int, required=True, help="Maximum number of leaves to retain")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()
    prune_random_leaves(args.treefile, args.outfile, args.max_leaves, args.seed)
