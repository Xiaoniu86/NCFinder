# Features
- Detect bidirectional transcripts (regions where a gene is upstream/downstream of a cluster).
- Identify gene body overlaps (sense and antisense overlaps).
- Classify clusters into:
  - Bidirectional
  - Gene body overlaps
  - Antisense gene body overlaps
  - Intergenic regions
  - Canonical assignments
- Outputs results as CSV files for further analysis.

# Project Structure
``` bash
NCFinder/
├── src/
│   ├── interval_utils.py         # Utilities for creating interval trees
│   ├── bidirectional_detect.py   # Functions for bidirectional transcript detection
│   ├── overlap_detect.py         # Functions for gene body overlap detection
│   ├── main.py                   # Main script to run the pipeline
│   └── __init__.py               # Package initialization
│
├── tests/
│   ├── test_bidirectional.py     # Unit tests for bidirectional detection
│   ├── test_overlap.py           # Unit tests for overlap detection
│   └── test_main.py              # Integration tests
│
├── data/
│   └── saccharomyces_cerevisiae_exons.csv   # Example input data
│
├── output/                       # Folder where result files are stored
│   └── (Generated CSV outputs)
│
├── .gitignore                    # Ignores unnecessary files
├── README.md                     # Documentation (this file)
├── requirements.txt              # Project dependencies
└── LICENSE                       # License information
