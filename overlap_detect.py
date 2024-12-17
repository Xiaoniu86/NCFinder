import pandas as pd
from intervaltree import IntervalTree

def overlap_detect_chromosome(cl_df, gene_df, chromosome):
    cl_data = cl_df[cl_df['chr'] == chromosome]
    gene_data = gene_df[gene_df['chromosome'] == chromosome]

    gene_tree = IntervalTree()
    for _, row in gene_data.iterrows():
        gene_tree[row['gene_start']:row['gene_end']] = (row['gene_name'], row['strand'])

    results = []
    for _, row in cl_data.iterrows():
        overlapping_intervals = gene_tree[row['start']:row['end']]
        results.append({
            **row,
            'overlap_gene': ",".join(set([i.data[0] for i in overlapping_intervals])),
            'antisense_gene': ",".join(set([i.data[0] for i in overlapping_intervals if i.data[1] != row['strand']]))
        })
    return pd.DataFrame(results)
