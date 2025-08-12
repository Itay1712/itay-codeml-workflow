import sys
from ete3 import Tree
import os
import re


## strip a trailing "#<digits>" from a leaf name
_HASH_SUFFIX = re.compile(r"#\d+$")
def strip_marking(name: str) -> str:
    return _HASH_SUFFIX.sub("", name)

def prune_tree(tree_file, taxa_to_remove):
    tree = Tree(tree_file, format=1)

    ## anything that is not found in the taxa to remove - keep
    taxa_to_keep = [
        leaf.name
        for leaf in tree.get_leaves()
        if strip_marking(leaf.name) not in taxa_to_remove
    ]
    if not taxa_to_keep:
        print("Warning: pruning removed all taxa; writing an empty/degenerate tree.")

    # all_taxa = set(tree.get_leaf_names())
    # print("num of taxa, before and after:", len(all_taxa), len(taxa_to_keep))
    # print("all taxa:", repr(list(all_taxa)))
    # print("to remove:", repr(taxa_to_remove))

    tree.prune(taxa_to_keep, preserve_branch_length=True)

    base_name, ext = os.path.splitext(tree_file)
    output_file = f"{base_name}_pruned{ext}"
    tree.write(outfile=output_file, format=1)
    print(f"Pruned tree saved as '{output_file}'")

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 prune_leaves_by_name.py tree_file taxon1 taxon2 ...")
        sys.exit(1)

    tree_file = sys.argv[1]
    taxa_to_remove = sys.argv[2:]
    prune_tree(tree_file, taxa_to_remove)

if __name__ == "__main__":
    main()
