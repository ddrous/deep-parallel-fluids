
## A simple rectangle regularly discretised
mutable struct Geometry
    Lx::Float64
    Ly::Float64
    Nx::Int64
    Ny::Int64
    b_ids::NamedTuple
end
geom = Geometry(1,1,8,5, (s=1, e=2, n=3, w=4))

# exit()

mutable struct Laplace
    bc:: Float64
    u:: Matrix{Float64}
    f:: Matrix{Float64}
end

## Go from vector of size (N)*(M) to appropriate matrix 
vec2mat(U::Vector{T}, N::Int64, M::Int64) where {T<:Number} = reverse(reshape(U, (N,M))', dims=1)

## Go from a matrix to a appropriate vector
mat2vec(u::Matrix{T}) where {T<:Number} = vec(permutedims(reverse(u, dims=1)))

function paddmat(u::Matrix{T}, bc::Vector{T,4}) where {T<:Number}
    ## Add 1 dimension and value of bc in each axis
    N, M = size(u)

    xvec1 = fill(bc[4], (1, M))
    xvec2 = fill(bc[2], (1, M))
    u = vcat(xvec1, u, xvec2)

    yvec1 = fill(bc[1], (N+2,1))
    yvec2 = fill(bc[3], (N+2,1))
    hcat(yvec2, u, yvec1)
end

function setbc(u::Matrix{T}, bc::Vector{T,4}, g:Geometry) where {T<:Number}
    ## Add 1 dimension and value of bc in each axis
    xvec = fill(bc[g.b_ids.w], (1, g.Ny))
    u = vcat(xvec,u)
    xvec = fill(bc[g.b_ids.e], (1, g.Ny))
    u = vcat(xvec,u)


    yvec = fill(bc[g.b_ids.s], (g.Nx+2,1))
    hcat(u, yvec)
    yvec = fill(bc[g.b_ids.n], (g.Nx+2,1))
    hcat(yvec, u)
end


function isonboundary(i::Int64, j::Int64, g::Geometry):Int64
    if j==1 
        return g.s
    elseif i==g.Nx+1
        return g.w
    elseif j==g.Ny+1
        return g.n
    elseif i==1
        return g.e
    else
        return 0
    end
end


function initlaplace(g::Geometry)
    bc = Vector{[0.0, 0.0, 0.0, 0.0]}
    u_small = zeros(Float64, (g.Nx-1,g.Ny-1))
    u = paddmat(u_small, bc)
    f = ones(Float64, (g.Nx-1,g.Ny-1))
    return Laplace(bc, u, f)
end

l = initlaplace(geom)

function solvelaplace(l::Laplace, g::Geometry)
    n = (g.Nx-1)*(g.Ny-1)
    U = zeros(Float64, n)
    F = vec2mat(l.f)
    A = zeros(Float64, (n,n))

    Δx = g.Lx/(g.Nx+1)
    Δy = g.Ly/(g.Ny+1)
    val1 = 1.0/(Δx^2)
    val2 = 1.0/(Δy^2)
    val3 = 2.0 * (val1+val2)

    bcvals = [val1, val2, val1, val2]

    for i in 2:g.Nx, j in 2:g.Ny
        i_glob = (i-1) + (j-2)*(N-1)
        A[i_glob,i_glob] = val3

        neighboor = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]
        for (ix, jy) in neighboor
            b_id = isonboundary(ix, jy, g.Nx, g.Ny)
            if b_id == 0
                ##Fill A
            else
                F[i_glob] += bcvals[b_id] * l.u[ix,jy]
            end
        end



    end


end