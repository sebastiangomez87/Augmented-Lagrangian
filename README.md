# Augmented-Lagrangian
Implementation of the augmented Lagrangian algorithm in Julia

For the equality problems
min f(x)
s.t. c(x)=0

the augmented Lagrangian, for a given multiplier u_i is given by (p is the penalization to deviations from the constraint):
L(x)=f(x)+u_i c(x)+p/2 ||c(x)||^2
The Lagrange multiplier is updated according to 
u_{i+1}=u_i+p c(x*)
where x* is the argmin that minimizes the augmented Lagrangian

For the inequality problems
min f(x)
s.t. c(x)>= 0

the augmented Lagrangian, for a given multiplier u_i is given by (p is the penalization to deviations from the constraint):
L(x)=f(x)+u_i c(x)+p/2 ||max{0,c(x)}||^2
The Lagrange multiplier is updated according to 
u_{i+1}=u_i+p max(0,c(x*))
where x* is the argmin that minimizes the augmented Lagrangian
