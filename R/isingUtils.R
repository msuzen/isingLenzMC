# 
#  Ising Model MC functions
#    Utility functions
#  (c) 2013 by Dr.Mehmet Suzen
#  GPLv3 or higher
# 

#
# Bexter 1982, eq 2.1.9
transferMatrix <- function(ikBt, J, H) {
  K  <- J*ikBt
  h  <- H*ikBt
  tm <- c(exp(K+h),exp(-K), exp(-K), exp(K-h))
  tm <- matrix(tm, 2, 2)
  list(tm=tm, evalues=eigen(tm)$values)
}
