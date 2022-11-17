## A simple rectangle regularly discretised
mutable struct Geometry
    Lx::Float64
    Ly::Float64
    Nx::Int64
    Ny::Int64
    Δx
    Δy
    x::Vector
    y::Vector
end



## A struct to hold our navierstokes problem attributes
mutable struct NavierStokes
    u:: Matrix{Float64}     ## Velocity
    p:: Matrix{Float64}     ## Pressure
    f:: Matrix{Float64}     ## External force
    ρ:: Float64             ## Density
end
