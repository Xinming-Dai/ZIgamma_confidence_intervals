# zigamma functions 
# ni0: the number of zero data points
# ni1: the number of non-zero data points
# zigdata: all data
# zigdata1: the section of non-zero data points

# convert a data frame contains zigamma data to a standard list format
zig = function(data){
  data0 = data[data == 0]
  data1 = data[data != 0]
  ni0 = length(data0)
  ni1 = length(data1)
  return(list(ni0 = ni0, ni1 = ni1, zigdata = data, zigdata1 = data1))
}

#-----------generate data from zero inflated gamma distribution---------------
rzigamma = function(n,delta,shape,scale){
  # n is the sample size.
  # delta is the proportion of 0 in binomial distribution.
  ni0=rbinom(1,n,prob=delta) # the prob of zeros which are generated using this way.
  ni1=n-ni0
  bdata=array(0,dim = ni0)
  gdata=rgamma(n=ni1,shape=shape,scale=scale)
  zigdata=c(bdata,gdata)
  return(list(ni0=ni0,ni1=ni1,zigdata=zigdata,zigdata1=gdata))
}