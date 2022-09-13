# Calculate Amino Acid Identity
Calculates the average amino acid identity among multiple taxa given their protein-coding genes. Application is currently directed towards (micro)eukaryotes and discerning deep phylogenomic relationships.

## Dependencies
[Python 3.6+](https://www.python.org/downloads/)\
[DIAMOND](https://github.com/bbuchfink/diamond)

## Usage
Requires a folder with translated protein sequences (one *per* taxon):
```
$ python3 aai_calc.py [FOLDER-WITH-AA-FASTAs]
```
Output includes a table summarizing the AAIs among all the taxa compared, as well as a substantially larger table with all pairwise comparisons used in the analyses. Note that this does solely incorporate Reciprocal Best Blast Hits (RBBHs).


## Planned Updates (V2)
- [ ] Update with phylogeny "awareness"
- [ ] Incorporate Ortholog assessment [OrthoFinder2](https://github.com/davidemms/OrthoFinder)
- [ ] Support for network analyses

### When using these tools in published research, please cite:
- Maurer-Alcalá XX, Gyaltshen Y, Heiss AA, Kim E. *In prep*. \"Sunset for Excavata: Analyses of Genome Architecture and Lineage-Specific Genes Break Apart the Contentious Supergroup"\.
