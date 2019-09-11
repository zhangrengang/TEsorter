# LTR_classifier
It is coded for [LTR_retriever](https://github.com/oushujun/LTR_retriever) to classify LTRs. It can also be used to classify any other TE sequences, including Class I and Class II elements which are covered by the [REXdb](http://repeatexplorer.org/?page_id=918) database.

### Installation ###
Dependencies:
+	[python 2.7](https://www.python.org/)
   +   [biopython](https://biopython.org/): quickly install by `pip install biopython`
   +   [parallel python](https://www.parallelpython.com/): quickly install by `pip install pp`
+	[hmmscan 3.1b2](http://hmmer.org/)
+   [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
 
`git clone https://github.com/zhangrengang/LTR_classifier`

### Quick Start ###
```
git clone https://github.com/zhangrengang/LTR_classifier
cd LTR_classifier
cd test

# run
python ../LTR_classifier.py rice6.9.5.liban
```
By default, the newly released [REXdb](http://repeatexplorer.org/?page_id=918) ([viridiplantae_v3.0 + metazoa_v3](https://bitbucket.org/petrnovak/re_databases)) database is used, which is more sensitive and more common and thus is recommended. 
[GyDB](http://gydb.org/) can also be used:
```
python ../LTR_classifier.py rice6.9.5.liban -db gydb
```
To speed up, use more processors [default=4]:
```
python ../LTR_classifier.py rice6.9.5.liban -p 36
```
To improve sensitivity, reduce the criteria:
```
python ../LTR_classifier.py rice6.9.5.liban -p 20 -cov 10 -eval 1e-2
```
To improve specificity, increase the criteria and disable the pass2 mode:
```
python ../LTR_classifier.py rice6.9.5.liban -p 20 -cov 30 -eval 1e-5 -dp2
```
To do with [TE polyprotein sequences](http://www.repeatmasker.org/RMDownload.html):
```
python ../LTR_classifier.py RepeatPeps.lib -st prot -p 20
```
### Outputs ###
```
rice6.9.5.liban.rexdb.domtbl        HMMScan raw output
rice6.9.5.liban.rexdb.faa           protein sequences of domain, which can be used for phylogenetic analysis.
rice6.9.5.liban.rexdb.dom.tsv       inner domains of LTRs, which might be used to filter domains based on their scores and coverages.
rice6.9.5.liban.rexdb.gff3          domain annotations
rice6.9.5.liban.rexdb.cls.tsv       TEs/LTRs classifications
	Column 5: "yes" means one LTR Copia/Gypsy element with full GAG-POL domains.
rice6.9.5.liban.rexdb.classified    library for RepeatMakser
```

### Usage ###
```
$ python LTR_classifier.py  -h
usage: LTR_classifier.py [-h] [-v]
                         [-db {rexdb,rexdb-plant,rexdb-metazoa,gydb}]
                         [-st {nucl,prot}] [-pre PREFIX] [-fw] [-p PROCESSORS]
                         [-tmp TMP_DIR] [-cov MIN_COVERAGE] [-eval MAX_EVALUE]
                         [-dp2] [-nolib] [--no-cleanup]
                         sequence

positional arguments:
  sequence              input TE sequences in fasta format [required]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -db {rexdb,rexdb-plant,rexdb-metazoa,gydb}, --hmm-database {rexdb,rexdb-plant,rexdb-metazoa,gydb}
                        the database used [default=rexdb]
  -st {nucl,prot}, --seq-type {nucl,prot}
                        'nucl' for DNA or 'prot' for protein [default=nucl]
  -pre PREFIX, --prefix PREFIX
                        output prefix [default='{-s}.{-db}']
  -fw, --force-write-hmmscan
                        if False, will use the existed hmmscan outfile and
                        skip hmmscan [default=False]
  -p PROCESSORS, --processors PROCESSORS
                        processors to use [default=4]
  -tmp TMP_DIR, --tmp-dir TMP_DIR
                        directory for temporary files [default=./tmp]
  -cov MIN_COVERAGE, --min-coverage MIN_COVERAGE
                        mininum coverage for protein domains in HMMScan
                        output.[default=20]
  -eval MAX_EVALUE, --max-evalue MAX_EVALUE
                        maxinum E-value for protein domains in HMMScan
                        output.[default=0.001]
  -dp2, --disable-pass2
                        do not further classify the unclassified sequences
                        [default=False for `nucl`, True for `prot`]
  -nolib, --no-library  do not generate a library file for RepeatMasker
                        [default=False]
  --no-cleanup          do not clean up the temporary directory
                        [default=False]
```

### Limitations ###
1. For each domain (e.g. RT), only the best hit with the highest score will output, which means: 1) if frame is shifted, only one part can be annotated; 2) for example, if two or more RT domains are present in one query sequence, only one of these RT domains will be annotated.
2. Many LTRs cannot be classified due to no hit, which might be because: 1) the database is still incompleted; 2) some LTRs may have too many mutations such as frame shifts and stop gains; 3) some LTRs may be false positive. For the test data set ([rice6.9.5.liban](https://raw.githubusercontent.com/oushujun/EDTA/master/database/rice6.9.5.liban)), ~84% LTRs (INT sequences) are classified.

### Further phylogenetic analyses ###
You may want to use the RT domains to analysis relationships among retrotransposons (LTR, LINE, DIRS, etc.). Here is an example:
```
# to extract RT domain sequences
cat rice6.9.5.liban.rexdb.dom.tsv | grep RT | python ../bin/get_record.py -i rice6.9.5.liban.rexdb.faa -o rice6.9.5.liban.rexdb.RT.faa -t fasta

# to align with MAFFT or other tools
mafft --auto rice6.9.5.liban.rexdb.RT.faa > rice6.9.5.liban.rexdb.RT.faa.aln

# to reconduct the phylogenetic tree with IQTREE or other tools
iqtree -s rice6.9.5.liban.rexdb.RT.faa.aln -bb 1000 -nt AUTO 

# Finally, visualize and edit the tree 'rice6.9.5.liban.rexdb.RT.faa.aln.treefile' with FigTree or other tools.
```
