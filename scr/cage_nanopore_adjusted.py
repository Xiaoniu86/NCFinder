import pandas as pd

###############################################
# 0) READ THE INPUT CSV
###############################################
df = pd.read_csv("./all/scer/decap/nanopore/tama/merged_gtf_dcp2_2_20_merge.csv")
df['Processed'] = False  # new column to mark used rows

###############################################
# STEP 1: 2-ROW PAIRS (CAGE + NANO) BY (ID, SEGMENT)
###############################################

def same_id_segment_cage_and_nanopore(group):
    """
    Return True if (within this ID+segment group) there are exactly 2 rows:
    one with 'dcp2cage' and one with 'dcp2nanopore'.
    """
    if len(group) != 2:
        return False
    return set(group['sourse']) == {'dcp2cage', 'dcp2nanopore'}

# Identify 2-row pairs
pairs_long_step1 = (
    df[~df['Processed']]  # only consider unprocessed
    .groupby(['ID', 'segment'], group_keys=False)
    .filter(same_id_segment_cage_and_nanopore)
)

def merge_cage_nano_pair(group):
    """
    Merges the 2-row group: (1) dcp2cage, (2) dcp2nanopore.
    """
    if len(group) != 2:
        return None
    
    cage = group[group['sourse'] == 'dcp2cage'].iloc[0]
    nano = group[group['sourse'] == 'dcp2nanopore'].iloc[0]

    strand = cage['stand']
    chromosome = cage['chr']
    ID_val = cage['ID']

    if strand == '+':
        new_start = cage['start']
        new_end   = nano['end']
    else:  # strand == '-'
        new_start = nano['start']
        new_end   = cage['end']

    new_length = abs(new_end - new_start) + 1
    new_segment = f"{cage['segment']};{nano['segment']}"
    new_annotation = f"{cage['annotation']};{nano['annotation']}"

    return pd.Series({
        'chr':        chromosome,
        'start':      new_start,
        'end':        new_end,
        'stand':      strand,
        'Length':     new_length,
        'ID':         ID_val,
        'segment':    new_segment,
        'sourse':     'cage_nano_adjust',
        'annotation': new_annotation
    })

# Merge and collect results for Step 1
new_df_step1 = (
    pairs_long_step1
    .groupby(['ID','segment'], group_keys=False)
    .apply(merge_cage_nano_pair)
    .dropna()
    .reset_index(drop=True)
)

# Mark those rows as processed in df
df.loc[pairs_long_step1.index, 'Processed'] = True


###############################################
# STEP 2: 2-ROW PAIRS (CAGE + NANO) WITH NUMERIC CAGE
###############################################

def two_segments_and_numeric_cage(g):
    """
    Return True if group g has exactly two rows:
      - One 'dcp2cage' with numeric annotation
      - One 'dcp2nanopore'
    """
    if len(g) != 2:
        return False
    if set(g['sourse']) != {'dcp2cage', 'dcp2nanopore'}:
        return False
    
    cage_row = g[g['sourse'] == 'dcp2cage'].iloc[0]
    return str(cage_row['annotation']).isdigit()

# Among unprocessed, identify such pairs
filtered_step2 = (
    df[~df['Processed']]
    .groupby('ID', group_keys=False)
    .filter(two_segments_and_numeric_cage)
)

# Merge them
new_rows_step2 = []
for ID_val, group in filtered_step2.groupby('ID'):
    cage_row = group[group['sourse'] == 'dcp2cage'].iloc[0]
    nano_row = group[group['sourse'] == 'dcp2nanopore'].iloc[0]
    
    new_chr = cage_row['chr']
    new_strand = cage_row['stand']
    if new_strand == '+':
        new_start = cage_row['start']
        new_end   = nano_row['end']
    else:
        new_start = cage_row['end']
        new_end   = nano_row['start']
    
    new_length = abs(new_end - new_start) + 1
    new_segment = f"{cage_row['segment']};{nano_row['segment']}"
    new_annotation = f"{cage_row['annotation']};{nano_row['annotation']}"

    new_rows_step2.append({
        'chr':        new_chr,
        'start':      new_start,
        'end':        new_end,
        'stand':      new_strand,
        'Length':     new_length,
        'ID':         ID_val,
        'segment':    new_segment,
        'sourse':     'cage_nano_adjust',
        'annotation': new_annotation
    })

new_df_step2 = pd.DataFrame(new_rows_step2, columns=[
    "chr","start","end","stand","Length","ID","segment","sourse","annotation"
])

# Mark these as processed
df.loc[filtered_step2.index, 'Processed'] = True


###############################################
# STEP 3: EXACT 3-ROW GROUPS
###############################################

# Load the tags file
df_tags = pd.read_csv("./all/scer/decap/annotation2024/merged_file_with_consensus0305.csv")
df_tags['annotation'] = df_tags['cluster'].astype('string')
df_tags = df_tags[['annotation','tags_dcp2']]

df_unprocessed_step3 = df[~df['Processed']].copy()

# Separate cage vs nano
cage_ = df_unprocessed_step3[df_unprocessed_step3['sourse'] == 'dcp2cage']
nano_ = df_unprocessed_step3[df_unprocessed_step3['sourse'] == 'dcp2nanopore']

# Merge cage rows with tags
cage_merged = pd.merge(cage_, df_tags, on='annotation', how='left').dropna(subset=['tags_dcp2'])
df_merged_step3 = pd.concat([cage_merged, nano_], ignore_index=True)

def pick_three_segments(g):
    """
    For a group of exactly 3 rows:
      - if 2 cage + 1 nano => pick cage w/ highest tags_dcp2
      - if 1 cage + 2 nano => pick best nano by coordinate logic
    """
    if len(g) != 3:
        return None
    
    cage_rows = g[g['sourse'] == 'dcp2cage']
    nano_rows = g[g['sourse'] == 'dcp2nanopore']

    # Filter cage rows to numeric annotation
    cage_rows = cage_rows[cage_rows['annotation'].str.isdigit()]

    if len(cage_rows) == 0 or len(nano_rows) == 0:
        return None
    
    strands = g['stand'].unique()
    if len(strands) != 1:
        return None
    strand = strands[0]

    chrs = g['chr'].unique()
    if len(chrs) != 1:
        return None
    chromosome = chrs[0]

    def build_new_row(cage, nano):
        if strand == '+':
            ns = cage['start']
            ne = nano['end']
        else:
            ns = cage['end']
            ne = nano['start']

        length_ = abs(ne - ns) + 1
        seg_ = f"{cage['segment']};{nano['segment']}"
        ann_ = f"{cage['annotation']};{nano['annotation']}"
        tags_ = cage.get('tags_dcp2', None)
        return pd.Series({
            'chr':        chromosome,
            'start':      ns,
            'end':        ne,
            'stand':      strand,
            'Length':     length_,
            'ID':         cage['ID'],
            'segment':    seg_,
            'sourse':     'cage_nano_adjust',
            'annotation': ann_,
            'tags_dcp2':  tags_
        })

    # Case A: 2 cage + 1 nano
    if len(cage_rows) == 2 and len(nano_rows) == 1:
        chosen_cage = cage_rows.loc[cage_rows['tags_dcp2'].idxmax()]
        chosen_nano = nano_rows.iloc[0]
        return build_new_row(chosen_cage, chosen_nano)

    # Case B: 1 cage + 2 nano
    if len(cage_rows) == 1 and len(nano_rows) == 2:
        cage_row = cage_rows.iloc[0]
        cage_start = cage_row['start']
        cage_end   = cage_row['end']
        
        if strand == '+':
            bigger = nano_rows[nano_rows['start'] > cage_start]
            if len(bigger) == 1:
                chosen_nano = bigger.iloc[0]
            elif len(bigger) > 1:
                diffs = (bigger['start'] - cage_start).abs()
                chosen_nano = bigger.loc[diffs.idxmin()]
            else:
                # none bigger => pick the one w/ min end
                min_end_idx = nano_rows['end'].idxmin()
                chosen_nano = nano_rows.loc[min_end_idx]
        else:
            smaller = nano_rows[nano_rows['end'] < cage_end]
            if len(smaller) == 1:
                chosen_nano = smaller.iloc[0]
            elif len(smaller) > 1:
                diffs = (smaller['end'] - cage_end).abs()
                chosen_nano = smaller.loc[diffs.idxmin()]
            else:
                bigger_ = nano_rows[nano_rows['end'] > cage_end]
                if len(bigger_) == 0:
                    return None
                chosen_nano = bigger_.loc[bigger_['start'].idxmax()]
        
        return build_new_row(cage_row, chosen_nano)

    return None

result_step3 = (
    df_merged_step3
    .groupby('ID', group_keys=False)
    .apply(pick_three_segments)
    .dropna()
    .reset_index(drop=True)
)

# Mark the entire ID if used
ids_used_step3 = result_step3['ID'].unique()
df_to_mark_3 = df_unprocessed_step3[df_unprocessed_step3['ID'].isin(ids_used_step3)]
df.loc[df_to_mark_3.index, 'Processed'] = True

###############################################
# COMBINE STEPS 1-3 MERGED ROWS
###############################################
all_merged_rows = pd.concat([new_df_step1, new_df_step2, result_step3], ignore_index=True)

###############################################
# CREATE df_final (UPDATED) + merged_all
###############################################
df_final = df.copy()
merged_all = all_merged_rows.copy()

###############################################
# UNIQUE LABEL STEP (CAGE-ONLY OR NANO-ONLY GROUPS)
###############################################

def assign_unique_label(group):
    """
    Mark groups as 'cage' or 'nanopore' if they contain only that source type
    (in any variants), e.g. multiple 'dcp2nanopore' rows.
    Otherwise leave 'unique' blank and do not change 'Processed'.
    """
    source_set = set(group['sourse'])

    # Case 1: Entire group has no cage => all nanopore-like sources
    if 'dcp2cage' not in source_set and 'dcp2nanopore' in source_set:
        group['unique'] = 'nanopore'
        group['Processed'] = True

    # Case 2: Entire group has no nanopore => all cage-like sources
    elif 'dcp2nanopore' not in source_set and 'dcp2cage' in source_set:
        group['unique'] = 'cage'
        group['Processed'] = True

    else:
        group['unique'] = ''  # or None
        # do not change Processed in this else branch
    return group

df_final = df_final.groupby('ID', group_keys=False).apply(assign_unique_label)

###############################################
# STEP 4: MULTI-ROW (>=4) GROUPS
# PICK BEST PAIRS ITERATIVELY
###############################################

df_tags2 = pd.read_csv("./all/scer/decap/annotation2024/merged_file_with_consensus0305.csv")
df_tags2['annotation'] = df_tags2['cluster'].astype(str)
df_tags2 = df_tags2[['annotation','tags_dcp2']]

# Focus only on unprocessed & not 'unique' in {cage, nanopore}
df_unprocessed = df_final[
    (~df_final['Processed']) &
    (~df_final['unique'].isin(['cage','nanopore']))
].copy()

# Separate cage vs nano
cage_ = df_unprocessed[df_unprocessed['sourse'] == 'dcp2cage'].copy()
nano_ = df_unprocessed[df_unprocessed['sourse'] == 'dcp2nanopore'].copy()

# Merge tags into cage rows
cage_merged = pd.merge(cage_, df_tags2, on='annotation', how='left', suffixes=('','_tagfile'))
# cage_merged = cage_merged.dropna(subset=['tags_dcp2'])  # if you want to skip no-tag cages

df_unprocessed_merged = pd.concat([cage_merged, nano_], ignore_index=True)

# Ensure tags_dcp2 is numeric
if 'tags_dcp2' not in df_unprocessed_merged.columns:
    df_unprocessed_merged['tags_dcp2'] = 0
else:
    df_unprocessed_merged['tags_dcp2'] = df_unprocessed_merged['tags_dcp2'].fillna(0)

merged_records = []
used_indices = set()

def merge_cage_nano(cage_row, nano_row):
    strand = cage_row['stand']
    chromosome = cage_row['chr']
    if strand == '+':
        new_start = cage_row['start']
        new_end   = nano_row['end']
    else:
        new_start = nano_row['start']
        new_end   = cage_row['end']
    new_length = abs(new_end - new_start) + 1
    new_segment = f"{cage_row['segment']};{nano_row['segment']}"
    new_annotation = f"{cage_row['annotation']};{nano_row['annotation']}"
    return {
        'chr':        chromosome,
        'start':      new_start,
        'end':        new_end,
        'stand':      strand,
        'Length':     new_length,
        'ID':         cage_row['ID'],
        'segment':    new_segment,
        'sourse':     'cage_nano_adjust',
        'annotation': new_annotation,
        'tags_dcp2':  cage_row.get('tags_dcp2', 0)
    }

# Group by ID
for ID_val, group in df_unprocessed_merged.groupby('ID'):
    used_in_this_id = set()
    for strand_val, subgroup in group.groupby('stand'):
        cage_rows = subgroup[subgroup['sourse'] == 'dcp2cage'].copy()
        nano_rows = subgroup[subgroup['sourse'] == 'dcp2nanopore'].copy()
        # Sort cage by descending tags_dcp2
        cage_rows = cage_rows.sort_values(by='tags_dcp2', ascending=False)

        used_nano_indices = set()

        for cage_idx, cage_row in cage_rows.iterrows():
            if cage_idx in used_indices:
                continue
            available_nano = nano_rows[~nano_rows.index.isin(used_nano_indices)]
            if available_nano.empty:
                continue
            if strand_val == '+':
                distances = (available_nano['start'] - cage_row['start']).abs()
            else:
                distances = (available_nano['end'] - cage_row['end']).abs()
            closest_nano_idx = distances.idxmin()
            nano_candidate = available_nano.loc[closest_nano_idx]
            merged_row = merge_cage_nano(cage_row, nano_candidate)
            merged_records.append(merged_row)

            used_indices.add(cage_idx)
            used_indices.add(closest_nano_idx)
            used_nano_indices.add(closest_nano_idx)
            used_in_this_id.add(cage_idx)
            used_in_this_id.add(closest_nano_idx)

    # If any used in this ID => mark entire ID as processed
    if used_in_this_id:
        df_final.loc[df_final['ID'] == ID_val, 'Processed'] = True

merged_df_strategy1 = pd.DataFrame(merged_records)

# Combine these new merges with our existing merges (all_merged_rows)
all_merged_rows = pd.concat([merged_all, merged_df_strategy1], ignore_index=True)

###############################################
# SAVE RESULTS
###############################################
all_merged_rows.to_csv("./all/scer/decap/nanopore/tama/dcp2_2_20_newCluster2.csv", index=False)
df_final.to_csv("./all/scer/decap/nanopore/tama/merged_gtf_dcp2_2_20_merge_allGroupProcessed2.csv", index=False)

###############################################
# CHECK PROCESSED STATS
###############################################
num_processed = df_final['Processed'].sum()  # True counts as 1, False as 0
print(f"Number of processed rows: {num_processed}")

# Check if any ID is still fully unprocessed
all_unprocessed_mask = ~(df_final.groupby("ID")["Processed"].any())
all_unprocessed_ids = all_unprocessed_mask[all_unprocessed_mask].index

if len(all_unprocessed_ids) == 0:
    print("No ID group is entirely unprocessed.")
else:
    print("The following IDs have all rows unprocessed:")
    print(all_unprocessed_ids.tolist())

print("\nDone! Full pipeline completed.")
