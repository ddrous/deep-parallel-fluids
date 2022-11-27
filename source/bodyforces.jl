include("geometry.jl")


## Application of body forces g using forward Euler
function forwardeuler(u::Array{Float64, 2}, v::Array{Float64, 2}, f::Array{Float64, 3}, Δt::Float64, g::Geometry)
    u_ = u .+ Δt .* f[:, 1:g.Ny, 1]
    v_ = v .+ Δt .* f[1:g.Nx, :, 2]
    return u_, v_
end
