using ForwardDiff, LinearAlgebra
function mynewton(f, x0; grad=x->ForwardDiff.gradient(f,x),hess=x->ForwardDiff.hessian(f,x),tol = 1e-6, maxiter = 1000)
    #This function uses the Newton method to solve min f(x)   
    x_old = x0
    normdiff = Inf
    iter = 1
    while normdiff > tol && iter <= maxiter
        #@show hess(x_old)
        x_new = x_old .- inv(hess(x_old))*grad(x_old)
        #@show x_new
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter = iter + 1
    end
    #return (value = x_old, funcRes=f(x_old), iter=iter)
    return x_old
end

function augL(f, g, h, x0; ρ0=7.0,Δρ=0.3, tol = 1e-6, maxiter = 1000)
    #This function uses aumented Lagrangian to solve min f(x) s.t. g(x)=0, h(x)>=0 where g(x) and h(x) are functions of x   
    λ_old=zeros(eltype(x0),length(g(x0))) #Multiplier for equality constraints
    θ_old=zeros(eltype(x0),length(h(x0))) #Multiplier for inequality constraints
    x_old=x0
    ρ_old=ρ0
    normdiff1 = Inf
    normdiff2 = Inf
    normdiff3 = Inf
    normdiff4 = Inf
    iter = 1
    while normdiff1 > tol && normdiff2 > tol && normdiff3 > tol && normdiff4 > tol && iter <= maxiter
        #augmentedLag(x)=f(x)+λ_old'*eqconst(x)-θ_old'*ineqconst(x)+ρ_old*(eqconst(x)'*eqconst(x)-ones(eltype(x0),length(ineqconst(x0)))'*ineqconst(x))
        z(x)=augmentedLag(f,g,h,λ_old,θ_old,ρ_old,x)
        #@show z(x_old)
        #sol=mynewton(z, x_old)
        #@show sol
        x_new=mynewton(z, x_old)
        #@show x_new
        normdiff1 = norm(x_new - x_old)
        x_old=x_new
        #@show g(x_new)
        #@show λ_old
        λ_new=λ_old.+ρ_old.*g(x_new)
        if minimum(θ_old.+ρ_old.*h(x_new))<0
            θ_new=zeros(eltype(x0),length(h(x0)))
        else
            θ_new=θ_old.+ρ_old.*h(x_new)
        end
        normdiff2 = norm(λ_new - λ_old)
        normdiff3 = norm(θ_new - θ_old)
        normdiff4 = norm(f(x_new) - f(x_old))
        ρ_old=ρ_old*(1+Δρ)
        iter=iter+1
    end
    return (value = x_old, funcRes=f(x_old), g=g(x_old),h=h(x_old), iter=iter)
end

function augmentedLag(f,g,h,λ_old,θ_old,ρ_old,x)
    if minimum(h(x))>0.0
        return f(x)+dot(λ_old,g(x))-dot(θ_old,h(x))+ρ_old*dot(g(x),g(x))
    else
        return f(x)+dot(λ_old,g(x))-dot(θ_old,h(x))+ρ_old*(dot(g(x),g(x))+dot(h(x),h(x)))
    end
    #return f(x)+dot(λ_old,g(x))+ρ_old*dot(g(x),g(x))
end