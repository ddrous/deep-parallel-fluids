include("geometry.jl")


## Application of body forces g using forward Euler
function forwardeuler(u::Array{Float64, 3}, f::Array{Float64, 3}, Δt::Float64, g::Geometry)
    u_ = u .+ Δt .* f
    return u_
end
