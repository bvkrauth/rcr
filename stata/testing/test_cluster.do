run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test CLUSTER option
*******************************************************************
quietly egen ID = seq()
quietly gen zero = 0
* No clustering
rcr SAT Small_Class ${controls}
savedresults save basic e()
* Clustering on group
rcr SAT Small_Class ${controls}, cluster(TCHID)
savedresults save cluster e()
assert reldif( e(betaxCI_H)  , 7.221434897597197 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 2.472053399759639 ) <  ${tol}
* Clustering on individual ID. Should give the same result as no clustering
rcr SAT Small_Class ${controls}, cluster(ID)
savedresults compare basic e(), tol(${tol})
* Clustering on a constant. Should give missing standard errors and
* (- inf, + inf) as the CI
rcr SAT Small_Class ${controls}, cluster(zero)
assert         e(betaxCI_L) < -8.99e305
assert         e(betaxCI_H) > -8.99e305
* Issue an error if the cluster variable does not exist
rcof "noisily rcr SAT Small_Class ${controls}, vce(cluster NOTHING) " == 111

* vce(cluster) is equivalent to cluster() and should produce the same result
rcr SAT Small_Class ${controls}, vce(cluster TCHID)
savedresults compare cluster e(), tol(${tol})
* Both cluster() and vce(cluster) can be provided if they are consistent
rcr SAT Small_Class ${controls}, vce(cluster TCHID) cluster(TCHID)
savedresults compare cluster e(), tol(${tol})
* Issue an error if the vce and cluster options are in conflict
rcof "noisily rcr SAT Small_Class ${controls}, vce(cluster TCHID) cluster(Girl)" == 100
