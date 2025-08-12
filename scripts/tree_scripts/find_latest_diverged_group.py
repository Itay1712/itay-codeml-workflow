#!/usr/bin/env python3
"""Determine the host group that diverged most recently in a rooted tree.

Given a Newick tree and a mapping file with at least ``Accession`` and ``Host``
columns, this script computes the distance from the root to every leaf in the
tree. For each host group, the maximal distance of its leaves from the root is
recorded. The host with the greatest such distance is considered the latest
diverged group and is printed to stdout.
"""

import sys
from typing import Dict


class Node:
    def __init__(self):
        self.children = []
        self.label = ""
        self.branch_length = None


def parse_newick(s: str) -> Node:
    """Parse a Newick string into a tree of ``Node`` objects."""

    def parse_subtree(i: int):
        node = Node()
        if s[i] == "(":
            i += 1
            while True:
                child, i = parse_subtree(i)
                node.children.append(child)
                if s[i] == ",":
                    i += 1
                elif s[i] == ")":
                    i += 1
                    break
                else:
                    raise ValueError(f"Unexpected character {s[i]!r} at position {i}")
            lbl, i = parse_label(i)
            node.label = lbl
        else:
            lbl, i = parse_label(i)
            node.label = lbl
        if i < len(s) and s[i] == ":":
            i += 1
            bl, i = parse_branch_length(i)
            node.branch_length = bl
        return node, i

    def parse_label(i: int):
        start = i
        while i < len(s) and s[i] not in [":", ",", ")", ";"]:
            i += 1
        return s[start:i], i

    def parse_branch_length(i: int):
        start = i
        while i < len(s) and s[i] not in [",", ")", ";"]:
            i += 1
        return s[start:i], i

    root, i = parse_subtree(0)
    if i < len(s) and s[i] == ";":
        i += 1
    return root


def build_depths(node: Node, current: float, depths: Dict[str, float]):
    label = node.label
    if label:
        base = label.split("|")[0]
        if base.endswith("#1"):
            base = base[:-2]
        depths[base] = current
    for child in node.children:
        bl = float(child.branch_length) if child.branch_length is not None else 0.0
        build_depths(child, current + bl, depths)


def find_latest_host(tree_file: str, mapping_file: str) -> str:
    with open(tree_file) as tf:
        tree_text = tf.read().strip()
    root = parse_newick(tree_text)
    depths: Dict[str, float] = {}
    build_depths(root, 0.0, depths)

    host_depth: Dict[str, float] = {}
    with open(mapping_file) as mf:
        header = mf.readline().strip().split("\t")
        try:
            acc_idx = header.index("Accession")
            host_idx = header.index("Host")
        except ValueError:
            raise SystemExit("Mapping file must contain 'Accession' and 'Host' columns")
        for line in mf:
            if not line.strip():
                continue
            parts = line.rstrip().split("\t")
            if len(parts) <= max(acc_idx, host_idx):
                continue
            acc = parts[acc_idx]
            host = parts[host_idx]
            depth = depths.get(acc)
            if depth is None:
                continue
            host_depth[host] = max(host_depth.get(host, float("-inf")), depth)
    if not host_depth:
        raise SystemExit("No host groups found in mapping file")
    return max(host_depth.items(), key=lambda x: x[1])[0]


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <tree_file> <mapping_file>")
    host = find_latest_host(sys.argv[1], sys.argv[2])
    print(host)
