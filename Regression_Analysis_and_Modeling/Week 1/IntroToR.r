############################################
##
##  A VERY Brief Introduction to R
##
##    Author: Ephraim M. Hanks
##            hanks@psu.edu 
##
############################################


####
####
#### (A) - Basic Operations in R  
####
####

##
## (A.1) Simple operations 
##

3+2
3-2
3*2
3/2
3^2
exp(3)
sqrt(3)
log(3)

##
## (A.2) Defining and storing scalars
##

a=4
b=3
a*b

##
## (A.3) Vectors
##

v=c(1,4,3,2)
v

u=1:4
u

## element-wise operations

v+2
v*2
u^2
u+v
u*v

##
## (A.4) Subsetting vectors
##

## subsetting
v
v[2]
v[1:3]
v[c(1,4)]

## logical arguments
v
v>2
v==2
v>=2

## subsetting using logical arguments
v
which(v>2)
v[which(v>2)]

####
####
#### (B) Matrices and matrix operations
####
####

##
## (B.1) Creating matrices
##


A=matrix(1:12,nrow=4,ncol=3)
A

b=4:1
b

c=matrix(1:4,nrow=4,ncol=1)
c

d=1:3
d

ones=matrix(1,nrow=4,ncol=1)
ones

I=diag(1,nrow=3,ncol=3)
I

Diag.mat=diag(1:3)
Diag.mat

##
## (B.2) element-wise operations
##  

A
2*A
A/2
A^2

A
d
A+d
A+1:4

##
## (B.3) Matrix multiplication
##

f=matrix(c(1,2,0),ncol=1)

A
f
A%*%f

A
Diag.mat
A%*%Diag.mat

##
## (B.4) subsetting matrices
##

A
A[2,]
A[,2]
A[2,2]
A[1:2,]
A[,-1]

##
## (B.5) inverting matrices
##

## create a square matrix
W=t(A)%*%A
W
## find inverse of W
W.inv=solve(W)
## check
W.inv%*%W
W%*%W.inv

