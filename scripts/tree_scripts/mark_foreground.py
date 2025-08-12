#!/usr/bin/env python3
import sys
import csv
from collections import Counter, defaultdict
from ete3 import Tree
from ete3.parser.newick import NewickError
import os

# ---------- Helpers ---------- #

def load_tree_best_effort(tree_path: str) -> Tree:
    if not os.path.exists(tree_path):
        sys.exit(f"Tree file not found: {tree_path}")

    attempts = [
        dict(format=1, quoted_node_names=True),
        dict(format=0, quoted_node_names=True),
        dict(format=8, quoted_node_names=True),
        dict(format=1, quoted_node_names=False),
        dict(format=0, quoted_node_names=False),
    ]
    last_err = None
    for kw in attempts:
        try:
            return Tree(tree_path, **kw)
        except NewickError as e:
            last_err = e
        except Exception as e:
            last_err = e
    sys.exit(
        "Failed to read tree. Last error:\n"
        f"{last_err}\n"
        "Tip: check Newick validity, node-name quoting, and try different 'format' flags."
    )

def read_mapping_tsv(path: str):
    rows = []
    with open(path, newline="") as fh:
        sniff = fh.readline()
        if not sniff:
            sys.exit(f"Empty mapping file: {path}")
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"Ancestral node", "Accession", "Host"}
        got = set(h.strip() for h in (reader.fieldnames or []))
        if not required.issubset(got):
            sys.exit(
                "Mapping header must contain exactly these columns:\n"
                "Ancestral node\tAccession\tHost"
            )
        for r in reader:
            anc = (r["Ancestral node"] or "").strip()
            acc = (r["Accession"] or "").strip()
            host = (r["Host"] or "").strip()
            if not (anc and acc and host):
                continue
            rows.append({"anc": anc, "acc": acc, "host": host})
    if not rows:
        sys.exit("No usable rows found in mapping file.")
    return rows

def index_by_name(t: Tree):
    idx = defaultdict(list)
    for n in t.traverse():
        nm = (n.name or "").strip()
        if nm:
            idx[nm].append(n)
    return idx

def resolve_node_by_name(t: Tree, name: str, name_index):
    if name.lower() == "root":
        return t.get_tree_root()
    nodes = name_index.get(name, [])
    if len(nodes) == 1:
        return nodes[0]
    elif len(nodes) > 1:
        print(f"[warn] Multiple nodes named '{name}'. Using the first encountered.", file=sys.stderr)
        return nodes[0]
    else:
        return None

def path_nodes_exclusive(ancestor, tip):
    path = []
    cur = tip
    while cur is not None and cur != ancestor:
        path.append(cur)
        cur = cur.up
    if cur != ancestor:
        return []
    path.reverse()
    return path

def majority_base_state_at_root(mapping_rows, root_name: str):
    root_like = []
    for r in mapping_rows:
        anc = r["anc"]
        if anc.lower() == "root" or anc == root_name:
            root_like.append(r["host"])
    if not root_like:
        return None
    return Counter(root_like).most_common(1)[0][0]

# ---------- Core annotation logic ---------- #

def annotate_states(tree_path: str, map_path: str):
    t = load_tree_best_effort(tree_path)
    mapping = read_mapping_tsv(map_path)
    name_idx = index_by_name(t)
    root = t.get_tree_root()

    for n in t.traverse():
        n.state = None
        n.prev_state = None
        n._is_change_point = False

    base = majority_base_state_at_root(mapping, (root.name or "").strip())
    if base is None:
        base = "Unknown"
        print("[warn] Could not infer a base state at root from mapping. Using 'Unknown'.", file=sys.stderr)
    root.state = base

    for r in mapping:
        anc = resolve_node_by_name(t, r["anc"], name_idx)
        if anc is None:
            print(f"[warn] Ancestral node '{r['anc']}' not found in tree. Skipping row {r}.", file=sys.stderr)
            continue
        acc_nodes = name_idx.get(r["acc"], [])
        if not acc_nodes:
            print(f"[warn] Accession (leaf) '{r['acc']}' not found in tree. Skipping row {r}.", file=sys.stderr)
            continue

        tip = next((c for c in acc_nodes if c.is_leaf()), acc_nodes[0])

        if not (tip == anc or tip in anc.get_descendants()):
            print(f"[warn] '{r['acc']}' is not under ancestor '{r['anc']}'. Skipping row {r}.", file=sys.stderr)
            continue

        anc._is_change_point = True

        segment = path_nodes_exclusive(anc, tip)
        if not segment:
            continue

        for node_on_path in segment:
            prev = getattr(node_on_path, "state", None)
            if prev is not None and prev != r["host"]:
                print(
                    f"[warn] Conflicting state at node '{node_on_path.name or '[unnamed]'}': "
                    f"was '{prev}', new '{r['host']}'. Overwriting.",
                    file=sys.stderr
                )
            node_on_path.state = r["host"]

    for n in t.traverse("preorder"):
        if n is root:
            continue
        if n.state is None:
            n.state = n.up.state

    for n in t.traverse("preorder"):
        n.prev_state = None if n is root else n.up.state

    for n in t.traverse():
        child_diff = any((c.state != n.state) for c in n.children)
        n.is_state_change = bool(n._is_change_point or child_diff)

    return t

def mark_target_state_nodes(t: Tree, target_state: str, suffix="#1", include_unnamed=False):
    for n in t.traverse():
        if getattr(n, "state", None) == target_state:
            nm = (n.name or "").strip()
            if nm:
                if not nm.endswith(suffix):
                    n.name = nm + suffix
            elif include_unnamed:
                ## add marking to unnamed nodes - does not add placeholders
                n.name = suffix

def print_report_to_stderr(t: Tree):
    unnamed_counter = 0
    def safe_name(n):
        nonlocal unnamed_counter
        if (n.name or "").strip():
            return n.name
        unnamed_counter += 1
        return f"__internal_{unnamed_counter}"

    for n in t.traverse("preorder"):
        nm = safe_name(n)
        state = n.state if n.state is not None else "Unknown"
        prev = n.prev_state if n.prev_state is not None else "None"
        is_leaf = "True" if n.is_leaf() else "False"
        is_change = "True" if getattr(n, "is_state_change", False) else "False"
        print(f"Node: {nm}", file=sys.stderr)
        print(f"  State: {state}", file=sys.stderr)
        print(f"  Previous state: {prev}", file=sys.stderr)
        print(f"  Is state-change: {is_change}", file=sys.stderr)
        print(f"  Is leaf: {is_leaf}\n", file=sys.stderr)

def write_tree_to_stdout(t: Tree):
    ## stdout of the Newick (redirected to the tree file)
    newick = t.write(format=1)
    print(newick)

# ---------- CLI ----------

def main(tree_path: str, map_path: str, target_state: str):
    t = annotate_states(tree_path, map_path)
    mark_target_state_nodes(t, target_state, suffix="#1", include_unnamed=True)
    ## human-readable report of each node to STDERR
    # print_report_to_stderr(t)
    write_tree_to_stdout(t)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: python annotate_states.py <tree.nwk> <mapping.tsv> <target_state>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
