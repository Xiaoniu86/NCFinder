from intervaltree import Interval, IntervalTree

def build_interval_tree_for_plus_genes(gene_data, upstream_gap, downstream_gap):
    tree = IntervalTree()
    for _, gene in gene_data.iterrows():
        start, end = gene['gene_start'] - upstream_gap, gene['gene_start'] + downstream_gap
        tree[start:end] = gene['gene_name']
    return tree

def build_interval_tree_for_minus_genes(gene_data, upstream_gap, downstream_gap):
    tree = IntervalTree()
    for _, gene in gene_data.iterrows():
        start, end = gene['gene_end'] - upstream_gap, gene['gene_end'] + downstream_gap
        tree[start:end] = gene['gene_name']
    return tree
