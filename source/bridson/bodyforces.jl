include("geometry.jl")


## Application of body forces g using forward Euler
function forwardeuler(u::Matrix, Δt::Float64, f::Float64, g::Geometry)
    u_ = u .+ Δt .* f
    return u_
end
