import os
import pandas as pd
from src.bidirectional_detect import bidirect_detect_chromosome
from src.overlap_detect import overlap_detect_chromosome

def main():
    assigned_folder = './all/scer/hinnebusch/assigned/'
    unassigned_folder = './all/scer/hinnebusch/unassigned/'
    output_folder = './output/'
    os.makedirs(output_folder, exist_ok=True)

    for unassigned_file in os.listdir(unassigned_folder):
        if unassigned_file.endswith(".unassignedClusters.txt"):
            cl_df = pd.read_csv(os.path.join(unassigned_folder, unassigned_file), sep="\t")
            gene_df = pd.read_csv('./data/saccharomyces_cerevisiae_exons.csv')

            # Bidirectional detection
            results_bi = bidirect_detect_chromosome(cl_df, gene_df, "chrI", 200, 50)

            # Overlap detection
            results_overlap = overlap_detect_chromosome(cl_df, gene_df, "chrI")

            # Save combined results
            results_bi.to_csv(f"{output_folder}/{unassigned_file}_bidirectional.csv", index=False)
            results_overlap.to_csv(f"{output_folder}/{unassigned_file}_overlap.csv", index=False)

if __name__ == "__main__":
    main()
