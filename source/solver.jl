include("structures.jl")
include("geometry.jl")
include("constraints.jl")
include("advection.jl")
include("bodyforces.jl")
include("projection.jl")

using IterativeSolvers

## Defines advection of quantity q over velocity field u
function advect(l::NavierStokes, Δt::Float64, g::Geometry)
    # q_ = finitedifference(q, u, Δt, g)
    # u_ = semilagrangian(l.u, l.u, gatherfield_x, Δt, g)
    # v_ = semilagrangian(l.v, l.v, gatherfield_y, Δt, g)
    u_ = semilagrangian_x(l.u, l.u, Δt, g)
    v_ = semilagrangian_y(l.v, l.v, Δt, g)
    l.u = u_
    l.v = v_
    return nothing
end


## Application of body forces g
function bodyforces(l::NavierStokes, Δt::Float64, g::Geometry)
    u_, v_ = forwardeuler(l.u, l.v, l.f, Δt, g)
    l.u = u_
    l.v = v_
    return nothing
end



## Calculates the right divergence-free u
# function project(u::Matrix, Δt::Float64)
function project(l::NavierStokes, Δt::Float64, g::Geometry)
    n = (g.Nx-2)*(g.Ny-2)
    # A = zeros(Float64, (n,n))
    A = spzeros(Float64, (n,n))
    F = zeros(Float64, (n))

    Nx, Ny = g.Nx, g.Ny
    Δx, Δy = g.Δx, g.Δy
    val1 = Δt/(l.ρ*(Δx^2))
    val2 = Δt/(l.ρ*(Δy^2))
    val3 = 2*(val1 + val2)

    for i in 2:g.Nx-1, j in 2:g.Ny-1        ## only for fluid points
        global_id = (i-1) + (j-2)*(g.Nx-2)
        A[global_id, global_id] = val3
        div_x = (l.u[i+1,j]-l.u[i,j]) / Δx
        div_y = (l.v[i,j+1]-l.v[i,j]) / Δy
        F[global_id] = -div_x - div_y ## negative divergence

        neighbors = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]
        for (ix, jy) in neighbors
            b_name = label(ix, jy, g)
            if b_name == "inside"   ##fluid
                if ix==i-1
                    A[global_id, global_id-1] = -val1
                elseif ix==i+1
                    A[global_id, global_id+1] = -val1
                elseif jy==j-1
                    A[global_id, global_id-(g.Nx-2)] = -val2
                elseif jy==j+1
                    A[global_id, global_id+(g.Nx-2)] = -val2
                end
            elseif b_name == "north" ## solid
                A[global_id, global_id] += -val1
                F[global_id] += -(1/Δx)*l.u[i+1,j]
            elseif b_name == "south" ## solid
                A[global_id, global_id] += -val1
                F[global_id] += +(1/Δx)*l.u[i,j]
            elseif b_name == "west" ## solid
                A[global_id, global_id] += -val2
                F[global_id] += +(1/Δy)*l.v[i,j+1]
            elseif b_name == "east" ## air/empty
                ## Do exactly nothing !! neither on A not the rhs
                # A[global_id, global_id] += -val2
                # F[global_id] += +(1/(2*Δy))*l.u[i,j+1,2]
            end
        end
    end

    ## Solve linear system
    # P = A\F
    P = cg(A, F)
    l.p[2:g.Nx-1, 2:g.Ny-1] = vectomat(P, g.Nx-2, g.Ny-2)


    # Print constructed quantities
    # show(stdout, "text/plain", A)
    # show(stdout, "text/plain", F)
    # show(stdout, "text/plain", P)
    # println()

    ## Expand pressure vector boundary condition to get all values for gradient update
    for i in 2:g.Nx, j in 2:g.Ny
        if (i==1)||(i==g.Nx) || (j==1)||(j==g.Ny)
            b_name = label(i, j, g)
            if b_name=="north"
                l.p[i,j] = l.p[i+1,j] + l.u[i+1,j]*(g.Δx * l.ρ)/Δt
            elseif b_name=="south"
                l.p[i,j] = l.p[i-1,j] - l.u[i,j]*(g.Δx * l.ρ)/Δt
            elseif b_name=="west"
                l.p[i,j] = l.p[i,j+1] + l.v[i,j+1]*(g.Δy * l.ρ)/Δt
            else b_name=="east"
                # l.p[i,j] = l.p[i,j-2] + l.u[i,j-1,2]*(2 * g.Δy * l.ρ)/Δt
                l.p[i,j] = 0
            end
        end
    end

    ## Update velocity with pressure gradient update
    u_next = similar(l.u)
    v_next = similar(l.v)
    for i in 2:g.Nx-1, j in 2:g.Ny-1
        ## easy
        grad_p_x = (l.p[i+1,j]-l.p[i,j])/(Δx)
        grad_p_y = (l.p[i,j+1]-l.p[i,j])/(Δy)
        u_next[i,j] = l.u[i,j] - (Δt/(l.ρ))*grad_p_x
        v_next[i,j] = l.v[i,j] - (Δt/(l.ρ))*grad_p_y
    end

    l.u = u_next
    l.v = v_next
    return nothing
end


## Solve a Laplave problem
function solvenavierstokes(l::NavierStokes, Δt::Float64, g::Geometry)

    # u = l.u0
    # for t in 0:Δt:l.T 
    # Δt = 5*min(l.Δx, l.Δx) / maximum(l.u)

    ## Advec
    advect(l, Δt, g)

    ## Body forces
    bodyforces(l, Δt, g)

    ## Project
    project(l, Δt, g)

    # return l.u, l.p
    return nothing
end
