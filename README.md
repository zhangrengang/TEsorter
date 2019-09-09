# LTR_classifier
It is coded for [LTR_retriever](https://github.com/oushujun/LTR_retriever) to classify LTRs. It can also be used to classify any other TE sequences, including Class I and Class II elements which are covered by the [REXdb](http://repeatexplorer.org/?page_id=918) database.

### Installation ###
Dependencies:
+	[python 2.7](https://www.python.org/)
   +   [biopython](https://biopython.org/)
+	[hmmscan 3.1b2](http://hmmer.org/)
 
`git clone https://github.com/zhangrengang/LTR_classifier`

### Quick Start ###
```
git clone https://github.com/zhangrengang/LTR_classifier
cd LTR_classifier
cd test

python ../LTR_classifier.py LTRlibAnn rice6.9.5.liban
python ../LTR_classifier.py Classifier rice6.9.5.liban.rexdb.gff3 > rice6.9.5.liban.rexdb.gff3.anno
```
By default, the newly released [REXdb](http://repeatexplorer.org/?page_id=918) ([viridiplantae_v3.0 + metazoa_v3](https://bitbucket.org/petrnovak/re_databases)) database is used, which is more sensitive and more common and thus is recommended. 
[GyDB](http://gydb.org/) can also be used for testing, like:
```
python ../LTR_classifier.py LTRlibAnn rice6.9.5.liban gydb
python ../LTR_classifier.py Classifier rice6.9.5.liban.gydb.gff3 gydb > rice6.9.5.liban.gydb.gff3.anno
```

### Outputs ###
```
rice6.9.5.liban.aa                  translated LTR sequences
rice6.9.5.liban.rexdb.domtbl        HMMScan raw output
rice6.9.5.liban.rexdb.tsv           inner domains of LTRs, which might be used to filter domains based on their scores and coverages.
rice6.9.5.liban.rexdb.faa           protein sequences of domain, which can be used for phylogenetics analysis.
rice6.9.5.liban.rexdb.gff3          domain annotations
rice6.9.5.liban.rexdb.gff3.anno     LTR classifications
	Column 5: "cmpl" means one LTR Copia/Gypsy element with full GAG-POL domains.
```

### Limitations ###
1. For each domain (e.g. RT), only the best hit with the highest score will output, which means: 1) if frame is shifted, only one part can be annotated; 2) for example, if two or more RT domains are present in one query sequence, only one of these RT domains will be annotated.
2. Many LTRs cannot be classified due to no hit, which might be because: 1) the database is still incompleted; 2) some LTRs may have too many mutations such as frame shifts and stop gains; 3) some LTRs may be false positive. For the test data set ([rice6.9.5.liban](https://raw.githubusercontent.com/oushujun/EDTA/master/database/rice6.9.5.liban)), ~84% LTRs (INT sequences) are classified.

