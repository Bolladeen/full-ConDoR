# Fast Constrained Dollo Solver

This repository contains a Python script to solve the Fast Constrained Dollo problem. The solution provides a phylogenetic tree in Newick format, CSV files with solutions, and DOT files for visualization.

## Input Files

You will need the following files as inputs:

1. **Mutation Matrix (`-i`)**: A CSV file with mutation data.
2. **Total Read Count Matrix (`-r`)**: CSV for the total read counts.
3. **Variant Read Count Matrix (`-v`)**: CSV for variant reads.
4. **Optional Files**:
    - **SNP List (`-s`)**: File containing SNPs.
    - **SNV List (`-s2`)**: File containing SNVs.
    - **Subclonal Mutations (`--subclonal_mutations`)**: YAML file specifying subclonal mutations.
    - **Copy Number Profiles (`--cnp`)**: CSV file for copy number profiles.
    - **Amplicon Metadata (`-m`)**: Metadata for mutation filtering.

---

## Running the Script

To run the script:

```bash
python3 fast-condor.py -i <mutation_matrix.csv> -r <total_reads.csv> -v <variant_reads.csv> -o <output_prefix>
```

Example with filtering and subclonal refinement:

```bash
python3 fast-condor.py -i input.csv -r total_reads.csv -v variant_reads.csv -m metadata.csv --scr -o results/output_prefix
```

## Parameters

- `-i`: Path to the mutation matrix CSV file (required).
- `-r`: Path to total read count matrix CSV (required).
- `-v`: Path to variant read count matrix CSV (required).
- `-s`: Path to SNP list file (optional).
- `-s2`: Path to SNV list file (optional).
- `-a`: False positive error rate (default: 0.03).
- `-b`: False negative error rate (default: 0.03).
- `--ado`: Allelic drop-out precision (default: 5).
- `-k`: Maximum allowed losses for SNVs (default: 2).
- `--pt`: Presence threshold for filtering (default: 0.85).
- `--vt0`: VAF threshold for homozygous mutations (default: 0.75).
- `--vt1`: VAF threshold for absent mutations (default: 0.25).
- `--trt`: Total reads threshold for reliable data (default: 10).
- `--mft`: Missing fraction threshold (default: 0.2).
- `--scr`: Enable subclonal refinement.
- `--subclonal_mutations`: YAML file of subclonal mutations (requires `--scr`).
- `--cnp`: Copy number profiles CSV (optional).

For a full list of parameters:
```bash
python3 fast-condor.py -h
```

## Output Files

- **`<output_prefix>_B.csv`**: Inferred mutation burden matrix.
- **`<output_prefix>_tree.dot`**: DOT format of the phylogenetic tree with cells.
- **`<output_prefix>_tree_without_cells.dot`**: DOT format of the phylogenetic tree excluding cells.
- **`<output_prefix>_tree.newick`**: Phylogenetic tree in Newick format.
