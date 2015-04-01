perl bp_search2gff.pl -i /rhome/cjinfeng/BigData/00.RD/GenomeAlign/Lastz/output/MSU7vsOBA/net_axt/Chr1.axt.chain.prenet.net.axt -o Chr1.gff -f axt --version 3 > log 2> log2 &
python Transposon_Age.py --maf Chr1.axt.chain.prenet.net.axt.maf
python test.py
python Transposon_Age.py --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --maf Chr1.axt.chain.prenet.net.axt.maf

echo "OGL, glaberrima"
python Transposon_Age.py --maf ../input/MSU7vsOGL_axt_maf --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.OGL.list > log 2> log2 &
python Transposon_Age.py --maf ../input/MSU7vsOGL_ortholog_maf --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.OGL.orth.list > log 2> log2 &
python Transposon_Age.py --maf ../input/MSU7vsOBA_ortholog_maf --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.OBA.orth.list > log 2> log2 &
python Transposon_Age.py --maf ../input/MSU7vsOPU_ortholog_maf --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.OPU.orth.list > log 2> log2 &
python Transposon_Age.py --maf ../input/MSU7vsOBR_ortholog_maf --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.OBR.orth.list > log 2> log2 &

echo "test for multiprocess"
python Transposon_Age.py --maf ../input/MSU7vsOGL_ortholog_maf --gff test10.gff --output test10.age.OGL.orth.list > log 2> log2 &
python Transposon_Age_mp.py --maf ../input/MSU7vsOGL_ortholog_maf --gff test10.gff --output test10.age.OGL.orth.mp.list > log 2> log2 &


echo "multiprocess run"
qsub age_mp.sh

