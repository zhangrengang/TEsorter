python ../LTR_classifier.py LTRlibAnn rice6.9.5.liban
python ../LTR_classifier.py Classifier rice6.9.5.liban.rexdb.gff3 > rice6.9.5.liban.rexdb.gff3.anno

python ../LTR_classifier.py LTRlibAnn rice6.9.5.liban gydb
python ../LTR_classifier.py Classifier rice6.9.5.liban.gydb.gff3 gydb > rice6.9.5.liban.gydb.gff3.anno

python ../LTR_classifier.py LTRlibAnn rice6.9.5.liban rexdb-plant
python ../LTR_classifier.py Classifier rice6.9.5.liban.rexdb-plant.gff3 rexdb-plant > rice6.9.5.liban.rexdb-plant.gff3.anno
python ../LTR_classifier.py LTRlibAnn rice6.9.5.liban rexdb-metazoa
python ../LTR_classifier.py Classifier rice6.9.5.liban.rexdb-metazoa.gff3 rexdb-metazoa > rice6.9.5.liban.rexdb-metazoa.gff3.anno

