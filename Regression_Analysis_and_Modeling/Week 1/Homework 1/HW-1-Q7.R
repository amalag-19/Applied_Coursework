## inputs
A=matrix(c(1,3,-1,9,5,-2,-1,-1,2,2,-3,1),4,3)
b=cbind(c(2,1,-1))

## part (a)
A%*%b

## part (b)
t(b)%*%b

## part(c)
b%*%t(b)

## part(d)
t(A)%*%A

## part(e)
sum(diag(t(A)%*%A))

## part(f)
solve(t(A)%*%A)





