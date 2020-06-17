%% A simple controllability check
load('../../Data/NSID'); %load the data
svd(ctrb(NSID.At,NSID.Bt))
rank(ctrb(NSID.At,NSID.Bt))
