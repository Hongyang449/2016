# add sca path; sca is required to load db
import sys
sys.path.append('/Users/hyangl/Documents/python/pySCA-master')
import scaTools as sca
import cPickle as pickle
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

db = pickle.load(open('/Users/hyangl/Documents/python/pySCA-master/Outputs/PF00503_full.db','rb'))
Dseq = db['sequence']  #the results of scaProcessMSA
Dsca = db['sca']       #the results of scaCore
Dsect = db['sector']   #the results of scaSectorID

# scaij; ro.r.matrix == matrix in R; or ro.r.['matrix']
cij = Dsca['Csca']
nr,nc = cij.shape
scaij_gt = ro.r.matrix(cij, nrow=nr, ncol=nc)

# residue number - according to 1TND
seq = ro.IntVector(Dseq['ats']) # try ro.r['as.integer'](Dseq['ats'])
ro.r.assign("scaij_gt", scaij_gt)
ro.r.assign("inds_gt", seq)
ro.r("save(scaij_gt, inds_gt, file='sca.RData', compress=FALSE)")


