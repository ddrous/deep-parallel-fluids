# using Base

## A simple rectangle regularly discretised
Base.@kwdef struct Geometry
    Lx::Float64
    Ly::Float64
    Nx::Int64
    Ny::Int64
    Δx::Float64 = Lx / Nx
    Δy::Float64 = Ly / Ny
    x::Vector{Float64} = collect(0:Δx:Lx)
    y::Vector{Float64} = collect(0:Δy:Ly)
end



## A struct to hold our navierstokes problem attributes
mutable struct NavierStokes
    u:: Array{Float64, 3}     ## Velocity
    p:: Matrix{Float64}     ## Pressure
    f:: Array{Float64, 3}     ## External force
    ρ:: Float64             ## Density
end
