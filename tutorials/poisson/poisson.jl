using Plots, SparseArrays 
plotlyjs()

## A simple rectangle regularly discretised
mutable struct Geometry
    Lx::Float64
    Ly::Float64
    Nx::Int64
    Ny::Int64
end


## Reshapes a vector into a Julia-ready matrix for computation
vectomat(U::Vector{T}, N::Int64, M::Int64) where {T<:Number} = reshape(U, (N,M))


## Flattens a matrix into a Julia-ready vector for computation
mattovec(u::Matrix{T}) where {T<:Number} = vec(u)


## Arranges a Julia-ready matrix into a 2D tensor suitable for human visualisation
matto2d(u::Matrix{T}) where {T<:Number} = reverse(permutedims(u), dims=1)


## A struct to hold out Poisson problem attributes
mutable struct Poisson
    bc:: Float64
    u:: Matrix{Float64}
    f:: Matrix{Float64}
end


## Sets a boundary condition all around a matrix
function setbc(u::Matrix{T}, bc::Float64, g::Geometry) where {T<:Number}
    u[1,:] .= bc
    u[g.Nx+1,:] .= bc
    u[:, 1] .= bc
    u[:, g.Ny+1] .= bc
    return nothing
end


## Initiate a Poisson problem
function initpoisson(g::Geometry; bc::Float64=1.0, fval::Float64=1.0)
    u = zeros(Float64, (g.Nx+1,g.Ny+1))
    setbc(u, bc, g)
    f = fill(fval, (g.Nx+1,g.Ny+1))
    return Poisson(bc, u, f)
end


## Deternimes location of a given point in a rectangle geometry
## N.B. Assumes matrix form (Ny,Nx), not 2D form (Nx,Ny) obtained by rotation
function isonboundary(i::Int64, j::Int64, g::Geometry):AbstractString
    if j==1 
        return "west"
    elseif i==g.Nx+1
        return "south"
    elseif j==g.Ny+1
        return "east"
    elseif i==1
        return "north"
    else
        return "inside"
    end
end


## Solve a Laplave problem
function solvepoisson(l::Poisson, g::Geometry)
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

    return nothing
end


## visualize a frame of the Poisson problem's solution
function visualizepoisson(l::Poisson, g::Geometry)
    ## Display solution in a 2D plane
    println("\n----Solution----")
    display(matto2d(l.u))
    println()

    ## Construct plot data: y axis first, then x axis (Julia convention)
    x=range(0,g.Lx,length=g.Nx+1)
    y=range(0,g.Ly,length=g.Ny+1)
    data = [l.u[i,j] for j in 1:g.Ny+1, i in 1:g.Nx+1]

    anns = [(0, 5, Plots.text("  N=$(g.Nx)")),
            (0, 4.7, Plots.text("M=$(g.Ny)")),
            (1, 5, Plots.text("Δx=$(g.Lx/g.Nx)")),
            (1, 4.7, Plots.text("Δy=$(g.Ly/g.Ny)")),
            (2, 5, Plots.text("f=$(l.f[1,1])")),
            (2, 4.7, Plots.text("   bc=$(l.bc)")),
            (0, 0, Plots.text("RDN")),
            (2, 5.6, Plots.text(" "))
            ]
    title = "Numerical solution of the Poisson problem"

    plot(x, y, data,
        xlabel="x", ylabel="y", zlabel="u", 
        camera=(30, 30), 
        st=:surface,
        annotations=anns,
        title=title,
        titlefontsize=18)
end


## Visiualize a gif of the Poisson problem's solution
function makepoissongif(c::Tuple, l::Poisson, g::Geometry)
    x=range(0, g.Lx, length=g.Nx+1)
    y=range(0, g.Ly, length=g.Ny+1)
    data = [l.u[i,j] for j in 1:g.Ny+1, i in 1:g.Nx+1]

    anns = [(0, 5, Plots.text("  N=$(g.Nx)")),
            (0, 4.7, Plots.text("M=$(g.Ny)")),
            (1, 5, Plots.text("Δx=$(g.Lx/g.Nx)")),
            (1, 4.7, Plots.text("Δy=$(g.Ly/g.Ny)")),
            (2, 5, Plots.text("f=$(l.f[1,1])")),
            (2, 4.7, Plots.text("   bc=$(l.bc)")),
            (0, 0, Plots.text("RDN")),
            (2, 5.6, Plots.text(" "))
            ]
    title = "Numerical solution of the Poisson problem"

    @gif for cx in vcat(c[1]:-2:c[2], c[2]:2:c[1])
        plot(x, y, data,
        xlabel="x", ylabel="y", zlabel="u", 
        camera=(cx, 30), 
        st=:surface,
        annotations=anns,
        title=title,
        titlefontsize=18)
    end
end


### Run the code in a function: no global variables !
function main()
    geo = Geometry(1,1,100,50)
    pb = initpoisson(geo; bc=1.0, fval=2.0)
    @time solvepoisson(pb, geo)
    # visualizepoisson(pb, geo)
    makepoissongif((80, -20), pb, geo)
end


main()
