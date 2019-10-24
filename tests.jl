#Some tests for the Augmente Lagrangian method. g(x) are equality constraints while h(x) are inequaity constraints

function g(x)
    return 0
end

function h(x)
    return 0
end

function f(x)
    return (1.0-x[1])^2+1.0*(x[2]-x[1]^2)^2
end

#z(x)=augmentedLag(f,h,1.0,5.0,x)
#mynewton(z, [1.0; 1.0])
#mynewton(f,[1.0; 1.0])
@show augL(f, g, h, [-2.0, 3.0])


function g(x)
    return 0
end

function h(x)
    return x[1] + x[2] - 2.5
end

function f(x)
    return (1.0-x[1])^2+1.0*(x[2]-x[1]^2)^2
end

#z(x)=augmentedLag(f,h,1.0,5.0,x)
#mynewton(z, [1.0; 1.0])
#mynewton(f,[1.0; 1.0])
@show augL(f, g, h, [3.0, 0.0])
@show augL(f, g, h, [0.0, 0.0])

function g(x)
    return x[1]+x[2]-6.0
end

function h(x)
    return 0 #[x[1]; x[2]]
end

function f(x)
    return -(x[1]^0.5)*(x[2]^0.5)
end

@show augL(f, g, h, [1.0, 2.0])

function g(x)
    return x[1]-1
end

function h(x)
    return 0
end

function f(x)
    return x[1]^2+x[2]^2
end

@show augL(f, g, h, [2.0, 2.0],ρ0=700)


function g(x)
    return x[1]+x[2]-1
end

function h(x)
    return [x[1];x[2]]
end

function f(x)
    return x[1]^2+x[2]^2
end

@show augL(f, g, h, [0.0, 0.0])


function g(x)
    return x[1]+x[2]-1
end

function h(x)
    return [x[1];x[2]]
end

function f(x)
    return -x[1]-log(x[2])
end

@show augL(f, g, h, [3.0, 1.0])


function g(x)
    return x[1]+x[2]-10
end

function h(x)
    return [x[1];x[2]]
end

function f(x)
    return -x[1]-log(x[2])
end

@show augL(f, g, h, [3.0, 1.0])

function g(x)
    return -x[1]-x[2]+10
end

function h(x)
    return [x[1];x[2]]
end

function f(x)
    return -(x[1])*(x[2])
end

@show augL(f, g, h, [0.0, 1.0])


