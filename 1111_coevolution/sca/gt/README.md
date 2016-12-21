## sca results are caculated via pySCA; Here I convert the results from python into R

## 1.terminal: /Users/hyangl/Documents/python/pySCA-master/

# directly use PF00503_rd2.an provided by them; try to generate it by myself from pfam alignment later.
./annotate_MSA.py Inputs/PF00503_full.txt -o Outputs/PF00503_full.an -a 'pfam'
./scaProcessMSA.py Outputs/PF00503_full.an -s 1TND -c B -f 'Bos taurus' -t -n
./scaCore.py Outputs/PF00503_full.db
./scaSectorID.py Outputs/PF00503_full.db

reference sequence index is: 2779
F7B680_HORSE|4|Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Laurasiatheria; Perissodactyla; Equidae; Equus.|0
-----EKHSRELEKKLKEDAEKDARTVKLLLLGAGESGKSTIVKQMKIIHQDGYSLEECLEFIAIIYGNTLQSILAIVRAMTTLNIQY--GDSARQDDARKLMHMADTI---EEGTMPKEMSDIIQRLWKDSGIQACFDRASEYQLNDSAGYYLSDLERLVTPGYVPTEQDVLRSRVKTTGIIETQFSFKDLNFRMFDVGGQRSERKKWIHCFEGVTCIIFIAALSAYDMVLVEDDEVNRMHESLHLFNSICNHRYFATTSIVLFLNKKDVFSEKIKKAHLSICFPDYNG-PNTYEDAGNYIKVQFLELNMR--RDVKEIYSHMTCATDTQNVKFVFDAVTDII--
truncating to reference sequence...
Conducting sequence and position filtering: alignment size is 4127 seqs, 313 pos
ATS and distmat size - ATS: 313, distmat: 313 x 313
Keeping 3136 sequences of 4127 sequences (after filtering for gaps)
Keeping 3102 sequences of 3136 sequences (after filtering for seq similarity)
After filtering: alignment size is 3102 seqs, 839 effective seqs, 310 pos
Final alignment parameters:
Number of sequences: M = 1258
Number of effective sequences: M' = 849    '
Number of alignment positions: L = 310
Number of positions in the ats: 310
Number of structure positions mapped: 310
Size of the distance matrix: 310 x 310
Opening database file Outputs/PF00503_full

## 2.python
rpy2_sca.py



