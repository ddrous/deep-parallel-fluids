
using Plots

## A simple rectangle regularly discretised
mutable struct Geometry
    Lx::Float64
    Ly::Float64
    Nx::Int64
    Ny::Int64
end

# exit()

mutable struct Laplace
    bc:: Float64
    u:: Matrix{Float64}
    f:: Matrix{Float64}
end

## Go from vector of size (N)*(M) to appropriate matrix 
vec2mat(U::Vector{T}, N::Int64, M::Int64) where {T<:Number} = reverse(permutedims(reshape(U, (M,N))), dims=1)

## Go from a matrix to a appropriate vector
mat2vec(u::Matrix{T}) where {T<:Number} = vec(permutedims(reverse(u, dims=1)))

function setbc(u::Matrix{T}, bc::Float64, g::Geometry) where {T<:Number}
    u[1,:] .= bc
    u[g.Nx+1,:] .= bc
    u[:, 1] .= bc
    u[:, g.Ny+1] .= bc
    return nothing
end


function isonboundary(i::Int64, j::Int64, g::Geometry):AbstractString
    if j==1 
        return "south"
    elseif i==g.Nx+1
        return "west"
    elseif j==g.Ny+1
        return "north"
    elseif i==1
        return "east"
    else
        return "inside"
    end
end


function initlaplace(g::Geometry)
    bc = 0.0
    u = zeros(Float64, (g.Nx+1,g.Ny+1))
    setbc(u, bc, g)
    f = ones(Float64, (g.Nx+1,g.Ny+1))
    return Laplace(bc, u, f)
end


function solvelaplace(l::Laplace, g::Geometry)
    n = (g.Nx-1)*(g.Ny-1)
    U = zeros(Float64, n)
    F = mat2vec(l.f[2:g.Nx, 2:g.Ny])
    A = zeros(Float64, (n,n))

    Δx = g.Lx/(g.Nx+1)
    Δy = g.Ly/(g.Ny+1)
    val1 = 1.0/(Δx^2)
    val2 = 1.0/(Δy^2)
    val3 = 2.0 * (val1+val2)

    bcvals = Dict("south"=>val1, "east"=>val2, "north"=>val1, "west"=>val2)

    for i in 2:g.Nx, j in 2:g.Ny
        global_id = (i-1) + (j-2)*(g.Nx-1)
        A[global_id,global_id] = val3

        neighbors = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]
        for (ix, jy) in neighbors
            b_name = isonboundary(ix, jy, g)
            if b_name == "inside"
                if ix==i-1
                    A[global_id-1,global_id] = -val1
                elseif ix==i+1
                    A[global_id+1,global_id] = -val1
                elseif jy==j-1
                    A[global_id,global_id-g.Nx+1] = -val2
                elseif jy==j+1
                    A[global_id,global_id+g.Nx-1] = -val2
                end
            else
                F[global_id] += bcvals[b_name] * l.u[ix,jy]
            end
        end
    end

    ## Solve system
    U = A\F
    l.u[2:g.Nx, 2:g.Ny] = vec2mat(U, g.Nx-1, g.Ny-1)

    return nothing

end

function visualizelaplace(l::Laplace, g::Geometry)
    x=range(0,1,length=g.Nx+1)
    y=range(0,1,length=g.Ny+1)

    data = [l.u[i,j] for j in 1:g.Ny+1, i in 1:g.Nx+1]
    plot(x, y, data, st=:surface, camera=(-30,30))
end





### Run the code

function main()
    geo = Geometry(1,1,8,5)
    pb = initlaplace(geo)

    solvelaplace(pb, geo)
    visualizelaplace(pb,geo)
end

main()

# exit()