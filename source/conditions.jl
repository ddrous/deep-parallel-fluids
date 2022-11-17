include("geometry.jl")

## Normal component is zeros
function nostick(u::Matrix, g::Geometry)
    for i in [1, g.Nx+1], j in [1, g.Ny+1]
        normal = boundarynormal(i, j, g)
        # dot(u[i,j], normal) == 0      ## Use a macro to force this condition !!
    end
    return nothing
end


## Normal and tangential components are zeros
function noslip(u::Matrix, g::Geometry)
    for i in [1, g.Nx+1], j in [1, g.Ny+1]
        u[i,j] = 0.0
    end
    return nothing
end


## Sets a boundary condition all around a matrix
function setbc(u::Matrix{T}, g::Geometry) where {T<:Number}
    noslip(u, g)
    return nothing
end
