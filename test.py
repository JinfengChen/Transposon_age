from gff3 import Gff3
gff = Gff3('Chr1.gff')
for line in gff.lines:
    for i in line.keys():
        print i, line[i]
