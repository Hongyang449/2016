# add sca path; sca is required to load db
import sys
sys.path.append('/Users/hyangl/Documents/python/pySCA-master')
import scaTools as sca
import cPickle as pickle
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

db = pickle.load(open('/Users/hyangl/Documents/python/pySCA-master/Outputs/PF00009_full.db','rb'))
Dseq = db['sequence']  #the results of scaProcessMSA
Dsca = db['sca']       #the results of scaCore
Dsect = db['sector']   #the results of scaSectorID

# scaij; ro.r.matrix == matrix in R; or ro.r.['matrix']
cij = Dsca['Csca']
nr,nc = cij.shape
scaij_eftu = ro.r.matrix(cij, nrow=nr, ncol=nc)

# residue number - according to 1TTT
seq = ro.IntVector(Dseq['ats']) # try ro.r['as.integer'](Dseq['ats'])
ro.r.assign("scaij_eftu", scaij_eftu)
ro.r.assign("inds_eftu", seq)
ro.r("save(scaij_eftu, inds_eftu, file='sca.RData', compress=FALSE)")


