include("structures.jl")
include("geometry.jl")
include("constraints.jl")
include("advection.jl")
include("bodyforces.jl")
include("projection.jl")

using IterativeSolvers

## Defines advection of quantity q over velocity field u
function advect(q::Array{Float64, 3}, u::Array{Float64, 3}, Δt::Float64, g::Geometry)
    # q_ = finitedifference(q, u, Δt, g)
    q_ = semilagrangian(q, u, Δt, g)
    return q_
end


## Application of body forces g
function bodyforces(u::Array{Float64, 3}, f::Array{Float64, 3}, Δt::Float64, g::Geometry)
    u_ = forwardeuler(u, f, Δt, g)
    return u_
end



## Calculates the right divergence-free u
# function project(u::Matrix, Δt::Float64)
function project(l::NavierStokes, Δt::Float64, g::Geometry)
    n = (g.Nx-1)*(g.Ny-1)
    # A = zeros(Float64, (n,n))
    A = spzeros(Float64, (n,n))
    F = zeros(Float64, (n))

    Δx = g.Lx/g.Nx
    Δy = g.Ly/g.Ny
    val1 = Δt/(4*l.ρ*(Δx^2))
    val2 = Δt/(4*l.ρ*(Δy^2))
    val3 = 2*(val1 + val2)

    for i in 2:g.Nx, j in 2:g.Ny
        global_id = (i-1) + (j-2)*(g.Nx-1)
        A[global_id, global_id] = val3
        div_x = (l.u[i+1,j,1]-l.u[i-1,j,1])/(2*Δx)
        div_y = (l.u[i+1,j,2]-l.u[i-1,j,2])/(2*Δy)
        F[global_id] = -div_x - div_y ## negative divergence

        neighbors = [(i-2, j), (i+2, j), (i, j-2), (i, j+2)]
        for (ix, jy) in neighbors
            b_name = is_on_or_beyond_boundary(ix, jy, g)
            if b_name == "inside"
                if ix==i-2
                    A[global_id, global_id-2] = -val1
                elseif ix==i+2
                    A[global_id, global_id+2] = -val1
                elseif jy==j-2
                    A[global_id, global_id-2*(g.Nx-1)] = -val2
                elseif jy==j+2
                    A[global_id, global_id+2*(g.Nx-1)] = -val2
                end
            elseif b_name == "north" ## solid
                A[global_id, global_id] += -val1
                F[global_id] += +(1/(2*Δx))*l.u[i+1,j,1]
            elseif b_name == "south" ## solid
                A[global_id, global_id] += -val1
                F[global_id] += -(1/(2*Δx))*l.u[i-1,j,1]
            elseif b_name == "west" ## solid
                A[global_id, global_id] += -val2
                F[global_id] += -(1/(2*Δy))*l.u[i,j-1,2]
            elseif b_name == "east" ## air/empty
                ## Do exactly nothing !! neither on A not the rhs
                # A[global_id, global_id] += -val2
                # F[global_id] += +(1/(2*Δy))*l.u[i,j+1,2]
            end
        end
    end


    ## Print constructed quantities
    # show(stdout, "text/plain", A)
    # show(stdout, "text/plain", F)

    ## Solve linear system
    # P = A\F
    P = cg(A, F)
    l.p[2:g.Nx, 2:g.Ny] = vectomat(P, g.Nx-1, g.Ny-1)

    ## Expand pressure vector boundary condition to get all values for gradient update
    for i in 2:g.Nx, j in 2:g.Ny
        if (i==1)||(i==g.Nx+1) || (j==1)||(j==g.Ny+1)
            b_name = is_on_or_beyond_boundary(i, j, g)
            if b_name=="north"
                l.p[i,j] = l.p[i+2,j] - l.u[i+1,j,1]*(2 * g.Δx * l.ρ)/Δt
            elseif b_name=="south"
                l.p[i,j] = l.p[i-2,j] + l.u[i-1,j,1]*(2 * g.Δx * l.ρ)/Δt
            elseif b_name=="west"
                l.p[i,j] = l.p[i,j+2] - l.u[i,j+1,2]*(2 * g.Δy * l.ρ)/Δt
            else b_name=="east"
                # l.p[i,j] = l.p[i,j-2] + l.u[i,j-1,2]*(2 * g.Δy * l.ρ)/Δt
                l.p[i,j] = 0
            end
        end
    end

    ## Update velocity with pressure gradient update
    u_next = similar(l.u)
    for i in 2:g.Nx, j in 2:g.Ny
        ## easy
        grad_p_x = (l.p[i+1,j]-l.p[i-1,j])/(2Δx)
        grad_p_y = (l.p[i,j+1]-l.p[i,j-1])/(2Δy)
        u_next[i,j,:] = l.u[i,j,:] - (Δt/(l.ρ)).*[grad_p_x, grad_p_y]
    end

    l.u = u_next
    return u_next
end


## Solve a Laplave problem
function solvenavierstokes(l::NavierStokes, Δt::Float64, g::Geometry)

    # u = l.u0
    # for t in 0:Δt:l.T 
    # Δt = 5*min(l.Δx, l.Δx) / maximum(l.u)

    ## Advec
    u_ = advect(l.u, l.u, Δt, g)

    ## Body forces
    u_ = bodyforces(u_, l.f, Δt, g)

    ## Project
    l.u = u_
    u = project(l, Δt, g)
    # end

    # return l.u, l.p
    return nothing
end
