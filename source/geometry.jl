include("structures.jl")

## Reshapes a vector into a Julia-ready matrix for computation
vectomat(U::Vector{T}, N::Int64, M::Int64) where {T<:Number} = reshape(U, (N,M))


## Flattens a matrix into a Julia-ready vector for computation
mattovec(u::Matrix{T}) where {T<:Number} = vec(u)


## Arranges a Julia-ready matrix into a 2D tensor suitable for human visualisation
matto2d(u::Matrix{T}) where {T<:Number} = reverse(permutedims(u), dims=1)


## Deternimes location of a given point in a rectangle geometry
## N.B. Assumes matrix form (Ny,Nx), not 2D form (Nx,Ny) obtained by rotation
function label(i::Int64, j::Int64, g::Geometry):AbstractString
    if j<=1
        return "west"   #west "solid"
    elseif i>=g.Nx
        return "south"  #south "solid"
    elseif j>=g.Ny
        return "east"   #east "empty"
    elseif i<=1
        return "north"  #north "solid"
    else
        return "inside" #inside "fluid"
    end
end


## Deternimes location of a given point in a rectangle geometry
## N.B. Assumes matrix form (Ny,Nx), not 2D form (Nx,Ny) obtained by rotation
function boundarynormal(i::Int64, j::Int64, g::Geometry):Vector
    if j==1
        return [0, -1]
    elseif i==g.Nx+1
        return [1, 0]
    elseif j==g.Ny+1
        return [0, 1]
    elseif i==1
        return [-1, 0]
    end
end


function gatherfield_xy(i::Int64, j::Int64, u::Matrix, v::Matrix):Float64
    u = (u[i,j] + u[i+1,j] + v[i,j] + v[i,j+1]) ./ 4.0
end

function gatherfield_x(i::Int64, j::Int64, u::Matrix):Float64
    u = (u[i,j] + u[i+1,j]) ./ 2.0
end

function gatherfield_y(i::Int64, j::Int64, u::Matrix):Float64
    u = (u[i,j] + u[i,j+1]) ./ 2.0
end
