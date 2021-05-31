# ``cORFusi`` - Correction of ORFs utilizing short-read information

## Summary
``cORFusi`` (correction of ORFs utilizing short-read information) corrects hybrid *de novo* assembled genomes using additional short-read assemblies and corresponding annotations.
Updating gene sequences and thus correcting ORFs leads to the annotation of complete genes.
Additional information regarding gene products and IDs could be added by [``Prokka``](https://github.com/tseemann/prokka) during annotation.

As input, an assembly (.fasta), an annotation (.gff) and a length threshold for up- and downstream region is required.
The genes listed in the annotation file are searched in the assembly via [``BLAST``](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
These candidates are filtered by two criteria: 1. 90% < similarity < 100%, 2. the length of the ``BLAST`` match may deviate by a maximum of 20% from the query length, i.e. length of the CDS.
Subsequently, upstream and downstream regions are defined for each candidate using the predefined length threshold, which are then blasted against the assembly.
Only exact matches are allowed.
This procedure determines the position of the ORF in the assembly.
If necessary, the reverse complement of the new sequence is calculated, which corresponds to cases 3) and 4), where the downstream region matched before the upstream region.
Finally, the sequence can be updated and the corrected assembly and corresponding log file are saved.

![Simplified overview of the ORF correction workflow with cORFusi.](workflow.pdf)

## Requirements
|Program/Package|Version|Note|
|---------------|-------|------|
|[``Python``](https://www.python.org/)|3.6.13||
|BCBio-gff|0.6.6||
|Biopython|1.76||
|[``BLAST``](https://blast.ncbi.nlm.nih.gov/Blast.cgi)|2.9.0|May work with other versions as well|

It is recommended to clone this repository and use a conda environment.
```bash
git clone https://github.com/sandraTriebel/corfusi
cd corfusi

conda env create -f corfusi.yaml  OR  conda create -n corfusi python=3.6.13 blast=2.9.0 bcbio-gff=0.6.6 biopython=1.76
conda activate corfusi
```

## Usage
After cloning the repository and installing all dependencies, we can start the ORF correction:
```bash
python corfusi.py -f assembly -g annotation -t int
```

|Short|Long|Description|
|-----|----|-----------|
|Mandatory parameters|
|**-h**|**--help** |Show this help message and exit|
|**-f**|**--fasta**|Path to input assembly file (.fasta)|   
|**-g**|**--gff**|Path to input annotation file (.gff)|
|**-t**|**--threshold**|Length threshold (integer) determing up- & downstream region|
|Optional parameters|
|**-p**|**--prefix**|Prefix for output files|
|**-o**|**--outdir**|Path to output folder|