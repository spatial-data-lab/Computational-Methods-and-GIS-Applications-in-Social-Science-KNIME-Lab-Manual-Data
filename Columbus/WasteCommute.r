require(lpSolve)
require(dplyr)

nREZ<-length(vec_res)
nEMP<-length(vec_emp)

#Set up constraint signs and right-hand sides
row.signs<-rep("=",nREZ)
row.rhs<-vec_res
col.signs<-rep("<=",nEMP)
col.rhs<-vec_emp


#Run to measure wasteful commuting
lpresult<-lp.transport(costs,"min",row.signs,row.rhs,col.signs,col.rhs)
print( paste('Value of objective function at optimum is ',lpresult$objval))
print( paste('The average optimum commute time  is ',round((lpresult$objval/sum(vec_res)),3)))