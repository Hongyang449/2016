## sca results are caculated via pySCA; Here I convert the results from python into R

## 1.terminal: /Users/hyangl/Documents/python/pySCA-master/

# directly use PF00009_rd2.an provided by them; try to generate it by myself from pfam alignment later.
./annotate_MSA.py Inputs/PF00009_full.txt -o Outputs/PF00009_full.an -a 'pfam'
./scaProcessMSA.py Outputs/PF00009_full.an -s 1TTT -c A -f 'Thermus aquatitus' -t -n
./scaCore.py Outputs/PF00009_full.db
./scaSectorID.py Outputs/PF00009_full.db

reference sequence index is: 14244
EFTU1_THET8|1|Bacteria; Deinococcus-Thermus; Deinococci; Thermales; Thermaceae; Thermus.|0
---KPHVNVGTIGHVDHGKTTLTAALTYVAAAENPNVEVKDY-GDIDKAPEERARGITINTAHVEYE-----TAKRHYSHVDCPGHADYIKNMITGAAQMDGAILVVSAADGPMPQTREHILLARQVGVPYIVVFMNKVDMVDDPELLDLVEMEVRDL-----------------------------------------------------------------------LNQYEFPGD-------EVPVIRGSALLALKIWELLDAIDEYIPT---
truncating to reference sequence...
Conducting sequence and position filtering: alignment size is 30330 seqs, 186 pos
ATS and distmat size - ATS: 186, distmat: 186 x 186
Keeping 28166 sequences of 30330 sequences (after filtering for gaps)
Keeping 27597 sequences of 28166 sequences (after filtering for seq similarity)
After filtering: alignment size is 27597 seqs, 4590 effective seqs, 167 pos
Final alignment parameters:
Number of sequences: M = 6885
Number of effective sequences: M' = 4457    '
Number of alignment positions: L = 167
Number of positions in the ats: 167
Number of structure positions mapped: 167
Size of the distance matrix: 167 x 167

## 2.python
rpy2_sca.py




