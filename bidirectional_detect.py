import pandas as pd
from src.interval_utils import build_interval_tree_for_plus_genes, build_interval_tree_for_minus_genes

def bidirect_detect_chromosome(cl_df, gene_df, chromosome, upstream_gap, downstream_gap):
    results = []
    # Separate data
    cl_minus = cl_df[(cl_df['chr'] == chromosome) & (cl_df['strand'] == '-')]
    gene_plus = gene_df[(gene_df['chromosome'] == chromosome) & (gene_df['strand'] == '+')]
    gene_plus_tree = build_interval_tree_for_plus_genes(gene_plus, upstream_gap, downstream_gap)

    cl_plus = cl_df[(cl_df['chr'] == chromosome) & (cl_df['strand'] == '+')]
    gene_minus = gene_df[(gene_df['chromosome'] == chromosome) & (gene_df['strand'] == '-')]
    gene_minus_tree = build_interval_tree_for_minus_genes(gene_minus, upstream_gap, downstream_gap)

    for _, cluster in cl_minus.iterrows():
        for overlap in gene_plus_tree[cluster['end']]:
            results.append({**cluster, 'bidirect_gene': overlap.data, 'gap': cluster['end']-overlap.begin})

    for _, cluster in cl_plus.iterrows():
        for overlap in gene_minus_tree[cluster['start']]:
            results.append({**cluster, 'bidirect_gene': overlap.data, 'gap': overlap.end-cluster['start']})
    return pd.DataFrame(results)
