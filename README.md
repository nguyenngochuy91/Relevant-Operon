# RGB: **R**elevant **G**ene **B**lock
## Purpose

RGB is a tool to find orthologous gene block to a reference gene block in projaryotic genomes. Gene blocks are genes co-located on the chromosome. In many cases, gene blocks are conserved between bacterial species, sometimes as operons, when genes are co-transcribed. The conservation is rarely absolute: gene loss, gain, duplication, block splitting and block fusion are frequently observed. 

RGB accepts a set of species and a gene block in a reference species. It then finds all gene blocks, orhtologous to the reference gene blocks.

RGB provides 2 method to find relevant gene block, naive method and approximated method. Naive method tries to exhaustively search all the combination, and 
approximated using greedy method that has an approximation result. 
## Requirements
* [Wget](https://www.gnu.org/software/wget/) 
* [Conda](https://conda.io/miniconda.html) (package manager so we don't have to use sudo)
* [Python 3+](https://www.python.org/download/releases/3.0/)
* [Biopython 1.63+](http://biopython.org/wiki/Download)
* [Clustalw](http://www.clustal.org/clustal2/#Download)
* [Muscle Alignment](https://www.drive5.com/muscle/downloads.htm)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [ETE3](http://etetoolkit.org/download/) (python framework for tree)
* [Mathplotlib](https://matplotlib.org/3.1.1/users/installing.html)
* [PDA](http://www.cibiv.at/software/pda/#download) (optional if you want to debias your tree base on Phylogenetic Diversity)

## Installation
Users can either use github interface Download button or type the following command in command line (assumming `git` was installed):
```bash
git clone https://github.com/nguyenngochuy91/Relevant-Operon.git
```
Install Miniconda (you can either export the path everytime you use RGB, or add it to the .bashrc file). Before using
the following command line, users will need to install [Wget](https://www.gnu.org/software/wget/). 
```bash
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O Miniconda-latest-Linux-x86_64.sh
bash Miniconda-latest-Linux-x86_64.sh -b -p ~/anaconda_ete/
export PATH=~/anaconda_ete/bin:$PATH;
conda update conda
```
We are going to create an environment to run **RGB**, replace `myenv` with your preference of environment name. 
```bash
conda create --name myenv
```
When conda asks you to proceed, type `y`. This will create an environment `myenv`. To activate this environment, type `source activate myenv`. To deactivate this environment, type  `source deactivate`

Now, we are going to install other dependencies within environment `myenv`. Activate the environemnt:
```bash 
source activate myenv
```

Install Biopython and ete3 using conda (highly recommended install biopython with conda)
```bash
conda install -c bioconda biopython ete3
```
Install ete_toolchain for visualization
```bash
conda install -c etetoolkit ete_toolchain
```

Install BLAST, ClustalW, MUSCLE 
```bash
conda install -c bioconda blast clustalw muscle
```

Install mathplotlib
```basg
conda install -c conda-forge matplotlib
```

For PDA, check installation instructions on this website: [PDA](http://www.cibiv.at/software/pda/#download)

## Usage

The easiest way to run the project is to execute the script [RGB](https://github.com/nguyenngochuy91/Relevant-Operon/blob/master/relevantOperon.py), which is inside the directory [Relevant-Operon]. 

### Run on example datasets
The users can run this script on the example data sets provided in directory [E_Coli](https://github.com/nguyenngochuy91/Relevant-Operon/tree/master/E_Coli) and [B_Sub](https://github.com/nguyenngochuy91/Relevant-Operon/tree/master/B_Sub). The two following command lines will run [relevantOperon](https://github.com/nguyenngochuy91/Relevant-Operon/blob/master/relevantOperon.py) on our 2 directories. 
The final results (visualization files) are stored in folder name that user choose (*result_naive*, and *result_approx* in the example)
#### E_Coli

Using naive method
```bash
./relevantOperon.py -g E_Coli/genomes/ -b E_Coli/gene_block_names_and_genes.txt -r NC_000913 -f E_Coli/phylo_order.txt -o result_naive -a N
```

Using approximated method
```bash
./relevantOperon.py -g E_Coli/genomes/ -b E_Coli/gene_block_names_and_genes.txt -r NC_000913 -f E_Coli/phylo_order.txt -o result_approx -a Y
```

We can then execute the script [analyze](https://github.com/nguyenngochuy91/Relevant-Operon/blob/master/analyze.py) to compare the runtime, deletion, duplication and split count
of our 2 method. The result is store in folder *analysis*, it contains file: Time(log10)Plot.png (running time in log10 scale), 3 graph file for pairwise event count, 3 graph file for event count versus reference gene block.
```bash
./analyze.py -n result_naive/E_Coli/ -a result_approx/E_Coli/ -b ./E_Coli/gene_block_names_and_genes.txt
```
```
usage: analyze.py [-h] [--naiveInput NAIVEINPUT] [--approxInput APPROXINPUT]
                  [--geneBlock GENEBLOCK]

optional arguments:
  -h, --help            show this help message and exit
  --naiveInput NAIVEINPUT, -n NAIVEINPUT
                        naive directory (either time or result)
  --approxInput APPROXINPUT, -a APPROXINPUT
                        aprrox directory (either time or result)
  --geneBlock GENEBLOCK, -b GENEBLOCK
                        gene_block_names_and_genes.txt file

```




#### B_Sub

Using naive method
```bash
./relevantOperon.py -g B_Sub/genomes/ -b B_Sub/gene_block_names_and_genes.txt -r NC_000964 -f B_Sub/phylo_order.txt -o result_naive -a N
```

Using approximated method
```bash
./relevantOperon.py -g B_Sub/genomes/ -b B_Sub/gene_block_names_and_genes.txt -r NC_000964 -f B_Sub/phylo_order.txt -o result_approx -a N
```
### Run on users' specific datasets
If the users wants to run the program on their own datasets, then they have to provide the following inputs:
  1. Directory that stores all the genomes file to study in genbank format. A file name should be something like this `NC_011567.gbk` where `NC_011567` is the locus name.
  2. Gene block text file that stores gene blocks in a reference species (this reference has to be in the genomes directory). The gene block format is tab delimited. The first column is the gene block name, then followed by the genes' name. For example, here is the `gene_block_names_and_genes.txt` file from Escheria coli K-12 MG1655.
```bash
astCADBE	astA	astB	astC	astD	astE
atpIBEFHAGDC	atpI	atpH	atpC	atpB	atpA	atpG	atpF	atpE	atpD
caiTABCDE	caiA	caiE	caiD	caiC	caiB	caiT
casABCDE12	casE	casD	casA	casC	casB	cas1	cas2
chbBCARFG	chbG	chbF	chbC	chbB	chbA	chbR
``` 
   3. Run RGB, the output is stored in directory `result`.
  ```bash
  ./relevantOperon.py -g genomes_directory -b gene_block_names_and_genes.txt -r ref_accession -a Y -o result
  ```
  ```
usage: relevantOperon.py [-h] [--genomes_directory GENOMES_DIRECTORY]
                         [--gene_blocks GENE_BLOCKS] [--reference REFERENCE]
                         [--filter FILTER] [--output OUTPUT] [--approx APPROX]

optional arguments:
  -h, --help            show this help message and exit
  --genomes_directory GENOMES_DIRECTORY, -g GENOMES_DIRECTORY
                        The directory that store all the genomes file
                        (E_Coli/genomes)
  --gene_blocks GENE_BLOCKS, -b GENE_BLOCKS
                        The E_Coli/gene_block_names_and_genes.txt file, this
                        file stores the operon name and its set of genes
  --reference REFERENCE, -r REFERENCE
                        The ncbi accession number for the reference genome
                        (NC_000913 for E_Coli and NC_000964 for B_Sub)
  --filter FILTER, -f FILTER
                        The filter file for creating the tree
                        (E_Coli/phylo_order.txt for E_Coli or
                        B_Sub/phylo_order.txt for B-Sub)
  --output OUTPUT, -o OUTPUT
                        Output directory to store the result
  --approx APPROX, -a APPROX
                        Using approx method (Y,N)


  ```
   
Besides, the users can also provide a filter text file. This filter file specifies the species to be included in the reconstruction analysis. The reason is that there might be families of species that are over representative in our genomes directory. This will reduce phylogenetic diversity and cause bias in our ancestral reconstruction. Hence, it is recomended to run [PDA](http://www.cibiv.at/software/pda/#download) on generated tree before proceeding further steps in our analysis. In order to achieve this, the user can follow the following instructions:
   1. Generate a phylogenetic tree from the genomes directory
   ```bash
   ./create_newick_tree.py -G genomes_directory -o tree_directory -f NONE -r ref_accession
   ```
   ```
   usage: create_newick_tree.py [-h] [-G DIRECTORY] [-o DIRECTORY] [-f FILE]
                             [-m STRING] [-t FILE] [-r REF] [-q]

optional arguments:
  -h, --help            show this help message and exit
  -G DIRECTORY, --genbank_directory DIRECTORY
                        Folder containing all genbank files for use by the
                        program.
  -o DIRECTORY, --outfolder DIRECTORY
                        Directory where the results of this program will be
                        stored.
  -f FILE, --filter FILE
                        File restrictiong which accession numbers this script
                        will process. If no file is provided, filtering is not
                        performed.
  -r REF, --ref REF     The reference genome number, such as NC_000913 for E_Coli
  -q, --quiet           Suppresses most program text outputs.

   ```
   2. Download and install [PDA](http://www.cibiv.at/software/pda/#download). Debias the phylogenetic tree using `PDA` program:
   ```bash
   ./debias.py -i tree_directory/out_tree.nwk -o pda_result.txt -s num -r ref_accession
   ```
   ```
   usage: debias.py [-h] [-i INPUT_TREE] [-o PDA_OUT] [-s TREE_SIZE] [-r REF]


optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_TREE, --input_tree INPUT_TREE
                        Input tree that we want to debias
  -o PDA_OUT, --pda_out PDA_OUT
                        Output of pda to be store.
  -s TREE_SIZE, --tree_size TREE_SIZE
                        Reduce the size of the tree to this size
  -r REF, --ref REF     Force to include the following species, here I force
                        to include the reference species

   ```
   3. Run RGB,  the output is stored in directory `result`. 
  ```bash
  ./relevantOperon.py -g genomes_directory -b gene_block_names_and_genes.txt -r ref_accession -f phylo_order.txt -a Y -o result
  ```




## Examples

Here are two gene blocks that were generated through our program. 
1. Gene block paaABCDEFGHIJK:

This gene block codes for genes involved in the catabolism of phenylacetate and it is not conserved between the group of studied bacteria.

![paaABCDEFGHIJK](https://github.com/nguyenngochuy91/Relevant-Operon/blob/master/paaABCDEFGHIJK.png "Gene block paaABCDEFGHIJK")
2. Gene block atpIBEFHAGDC:

This gene block catalyzes the synthesis of ATP from ADP and inorganic phosphate and it is very conserved between the group of studied bacteria.

![atpIBEFHAGDC](https://github.com/nguyenngochuy91/Relevant-Operon/blob/master/atpIBEFHAGDC.png "Gene block atpIBEFHAGDC")



## Credits
1. [An event-driven approach for studying gene block evolution in bacteria](http://bioinformatics.oxfordjournals.org/content/early/2015/04/13/bioinformatics.btv128.full)
2. [Tracing the ancestry of operons in bacteria](https://academic.oup.com/bioinformatics/article/35/17/2998/5300000)

