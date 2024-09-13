## fastppm-data

This repository provides code for constructing simulations
under the perfect phylogeny factorization model. 

The main simulation script is `bin/simulate.py`, which constructs clonal
trees, usage matrices, frequency matrices, and read count matrices 
for a specified number of mutations, samples, coverage, and random seed. 
The script has the following usage instructions:
```bash
usage: simulate.py [-h] --mutations MUTATIONS --samples SAMPLES --coverage COVERAGE --output OUTPUT [--seed SEED]

Simulate a clonal matrix, usage matrix and read counts.

options:
  -h, --help            show this help message and exit
  --mutations MUTATIONS
                        Number of mutations.
  --samples SAMPLES     Number of sequenced samples.
  --coverage COVERAGE   Expected sequencing coverage.
  --output OUTPUT       Output prefix.
  --seed SEED           Random seed.
```

We have wrapped this script in a Nextflow pipeline to allow for easy construction
of $240$ simulated datasets with $n = 100, 500, 1000, 2500$ mutations, $m = 3$ samples,
$c = 30, 100, 1000$ coverage, and $s = 1, ..., 20$ random seeds. 
The pipepine can be run by executing:
```bash
nextflow main.nf
```
This will create a directory `simulations` containing the simulated datasets. Each
simulated dataset is stored in a directory `simulations/n[MUTATIONS]_s[SAMPLES]_c[COVERAGE]_r[SEED]` 
and contains six files `sim_clonal_matrix.txt`, `sim_total_matrix.txt`, `sim_usage_matrix.txt`,
`sim_frequency_matrix.txt`, `sim_tree.txt`, and `sim_variant_matrix.txt`.
