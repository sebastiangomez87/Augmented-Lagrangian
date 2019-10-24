# Augmented-Lagrangian
Implementation of the augmented Lagrangian algorithm in Julia

For the equality problems
$min f(x)
s.t. c(x)=0$

the augmented Lagrangian, for a given multiplier $\lambda_i$ is given by ($\mu$ is the penalization to deviations from the constraint):
$L(x)=f(x)+\lambda_i c(x)+\frac{\mu}{2} ||c(x)||^2$
The Lagrange multiplier is updated according to 
$\lambda_{i+1}=\lambda_i+\mu c(x^{*})$
where $x^{*}$ is the argmin that minimizes the augmented Lagrangian
