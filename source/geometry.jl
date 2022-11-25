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
        return "solid"   #west
    elseif i>=g.Nx+1
        return "solid"  #south
    elseif j>=g.Ny+1
        return "empty"   #east
    elseif i<=1
        return "solid"  #north
    else
        return "fluid" #inside
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
