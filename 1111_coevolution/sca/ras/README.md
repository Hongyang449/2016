##########
## pfam ##
##########

## sca results are caculated via pySCA; Here I convert the results from python into R

## 1.terminal: /Users/hyangl/Documents/python/pySCA-master/

# directly use PF00071_rd2.an provided by them; try to generate it by myself from pfam alignment later.
./scaProcessMSA.py Inputs/PF00071_rd2.an -s 5P21 -c A -f 'Homo sapiens' -t -n
./scaCore.py Outputs/PF00071_rd2.db
./scaSectorID.py Outputs/PF00071_rd2.db

## 2.python
rpy2_sca.py

###########
## dimer ##
###########

## pySCA for ras dimer

# I use the ras subfamily MSA for the ras dimer predictions
./scaProcessMSA.py Inputs/ras_filtered.fa -s 5P21 -c A -f 'Homo sapiens'
./scaCore.py Outputs/ras_filtered.db
./scaSectorID.py Outputs/ras_filtered.db

sp|P01112|RASH_HUMAN/1-166
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQ
Len refseq 166, len refpos 166, Len alg seq 165, len pairalg 166, len gloalg 165
Conducting sequence and position filtering: alignment size is 1521 seqs, 165 pos
ATS and distmat size - ATS: 165, distmat: 165 x 165
Keeping 1521 sequences of 1521 sequences (after filtering for gaps)
Keeping 1521 sequences of 1521 sequences (after filtering for seq similarity)
After filtering: alignment size is 1521 seqs, 58 effective seqs, 162 pos
Final alignment parameters:
Number of sequences: M = 1521
Number of effective sequences: M' = 57'
Number of alignment positions: L = 162
Number of positions in the ats: 162
Number of structure positions mapped: 162
Size of the distance matrix: 162 x 162
Opening database file Outputs/ras_filtered





