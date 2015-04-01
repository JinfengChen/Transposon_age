#!/bin/bash
#PBS -l nodes=3:ppn=4
#PBS -l mem=12gb
#PBS -l walltime=100:00:00
#PBS -V

cd $PBS_O_WORKDIR

#python Transposon_Age_mp.py --maf ../input/MSU7vsOGL_ortholog_maf_clean --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.OGL.orth.clean.list > transposon.age.OGL.orth.clean.log 2> transposon.age.OGL.orth.clean.log2
#python Transposon_Age_mp.py --maf ../input/MSU7vsOID_ortholog_maf_clean --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.OID.orth.clean.list > transposon.age.OID.orth.clean.log 2> transposon.age.OID.orth.clean.log2
#python Transposon_Age_mp.py --maf ../input/MSU7vsONI_ortholog_maf_clean --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.ONI.orth.clean.list > transposon.age.ONI.orth.clean.log 2> transposon.age.ONI.orth.clean.log2
#python Transposon_Age_mp.py --maf ../input/MSU7vsORU_ortholog_maf_clean --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.ORU.orth.clean.list > transposon.age.ORU.orth.clean.log 2> transposon.age.ORU.orth.clean.log2
#python Transposon_Age_mp.py --maf ../input/MSU7vsOBA_ortholog_maf_clean --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.OBA.orth.clean.list > transposon.age.OBA.orth.clean.log 2> transposon.age.OBA.orth.clean.log2
python Transposon_Age_mp.py --maf ../input/MSU7vsOPU_ortholog_maf_clean --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.OPU.orth.clean.list > transposon.age.OPU.orth.clean.log 2> transposon.age.OPU.orth.clean.log2
python Transposon_Age_mp.py --maf ../input/MSU7vsOBR_ortholog_maf_clean --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.OBR.orth.clean.list > transposon.age.OBR.orth.clean.log 2> transposon.age.OBR.orth.clean.log2
python Transposon_Age_mp.py --maf ../input/MSU7vsHEG4_ortholog_maf_clean --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.HEG4.orth.clean.list > transposon.age.HEG4.orth.clean.log 2> transposon.age.HEG4.orth.clean.log2
python Transposon_Age_mp.py --maf ../input/MSU7vsA123_ortholog_maf_clean --gff ../input/MSU_r7.fa.RepeatMasker.out.gff --output transposon.age.A123.orth.clean.list > transposon.age.A123.orth.clean.log 2> transposon.age.A123.orth.clean.log2

echo "Done"
