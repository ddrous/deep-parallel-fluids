include("structures.jl")
include("geometry.jl")
include("conditions.jl")
include("advection.jl")
include("bodyforces.jl")
include("projection.jl")

## Defines advection of quantity q over velocity field u
function advect(q::Matrix, u::Matrix, Δt::Float64, g::Geometry)
    # q_ = finitedifference(q, u, Δt, g)
    q_ = semilagrangian(q, u, Δt, g)
    return q_
end


## Application of body forces g
function bodyforces(u::Matrix, Δt::Float64, f::Vector, g::Float64)
    u_ = forwardeuler(u, Δt, f, g)
    return u_
end

## Calculates the right divergence-free u
function project(u::Matrix, Δt::Float64)
    
    return u_next
end

## Solve a Laplave problem
function solvenavierstokes(l::NavierStokes, g::Geometry)
    n = (g.Nx-1)*(g.Ny-1)
    # A = zeros(Float64, (n,n))
    A = spzeros(Float64, (n,n))
    F = mattovec(l.f[2:g.Nx, 2:g.Ny])

    Δx = g.Lx/g.Nx
    Δy = g.Ly/g.Ny
    val1 = 1.0/(Δx^2)
    val2 = 1.0/(Δy^2)
    val3 = 2.0 * (val1+val2)

    bcvals = Dict("south"=>val1, "east"=>val2, "north"=>val1, "west"=>val2)

    for i in 2:g.Nx, j in 2:g.Ny
        global_id = (i-1) + (j-2)*(g.Nx-1)
        A[global_id, global_id] = val3

        neighbors = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]
        for (ix, jy) in neighbors
            b_name = isonboundary(ix, jy, g)
            if b_name == "inside"
                if ix==i-1
                    A[global_id, global_id-1] = -val1
                elseif ix==i+1
                    A[global_id, global_id+1] = -val1
                elseif jy==j-1
                    A[global_id, global_id-g.Nx+1] = -val2
                elseif jy==j+1
                    A[global_id, global_id+g.Nx-1] = -val2
                end
            else
                F[global_id] += bcvals[b_name] * l.u[ix,jy]
            end
        end
    end

    ## Print constructed quantities
    # show(stdout, "text/plain", A)
    # show(stdout, "text/plain", F)

    ## Solve linear system
    U = A\F
    l.u[2:g.Nx, 2:g.Ny] = vectomat(U, g.Nx-1, g.Ny-1)

    ## Advec

    ## Body forces

    ## Project

    return nothing
end
