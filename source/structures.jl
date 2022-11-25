# using Base

## A simple rectangle regularly discretised
Base.@kwdef struct Geometry
    Lx::Float64
    Ly::Float64
    Nx::Int64
    Ny::Int64
    Δx::Float64 = Lx / Nx
    Δy::Float64 = Ly / Ny
    x::Vector{Float64} = collect(1:1:Nx) .* Δx - (Δx / 2.0)
    y::Vector{Float64} = collect(1:1:Ny) .* Δy - (Δy / 2.0)
end



## A struct to hold our navierstokes problem attributes
mutable struct NavierStokes
    u:: Array{Float64, 2}     ## Velocity along x
    v:: Array{Float64, 2}     ## Velocity along y
    # w:: Array{Float64, 2}     ## Velocity along z
    p:: Array{Float64, 2}     ## Pressure
    f:: Array{Float64, 3}     ## External force
    ρ:: Float64             ## Density
end
