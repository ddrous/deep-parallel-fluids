include("geometry.jl")

## Numerical partial derivative
function numpartial(q_minus, q_mid, q_plus, space_step)
    # return (q_plus - q_mid) / space_step
    # return (q_mid - q_minus) / space_step
    return (q_plus - q_minus) / (2 * space_step)
end

## Defines advection of quantity q over velocity field u using finitedifference scheme
function finitedifference(q::Matrix, u::Matrix, Δt::Float64, g::Geometry)
    # (Nx_, Ny_) = size(q)
    q_ = similar(q)
    for i in 2:g.Nx, j in 2:g.Ny        ## Should be able to go from 1 to Nx+1
        partial_x = numpartial(q[i+1,j], q[i,j], q[i-1,j], g.Δx)
        partial_y = numpartial(q[i,j-1], q[i,j], q[i,j-1], g.Δy)
        q_[i,j] = q[i,j] - Δt * dot(u[i,j], [partial_x partial_y])
    end
    return q_
end



## Find the index left of xi
function locate(xi::Float64, x::Vector{Float64})
    n = size(x)
    i = length(x[ x .<= xi])
    if i == 0                ## If the point is on or before the left/down boundary
        return i, 0.0
    elseif i == n            ## If the point is on or after the right/up boundary
        return i-1, 1.0
    else
        α = (xi - x[i]) / (x[i+1] - x[i])
        return i, α
    end
end

## Interpolate between a and b, with barycenter coeff α
function interpolate(a, b, α)
    return (1 .- α) .* a + α .* b
end

## Defines advection of quantity q over velocity field u using semilagrangian scheme: 2nd order RK
function semilagrangian(q::Matrix, u::Matrix, Δt::Float64, g::Geometry)
    # (Nx_, Ny_) = size(q)
    q_ = similar(q)
    for i in 1:g.Nx+1, j in 1:g.Ny+1
        ## Find middle points
        x_mid, y_mid = [g.x[i]; g.y[j]] - 0.5*Δt*u[i,j]
        (i_mid, αi_mid) = locate(x_mid, g.x)        ## Could parallelise locate with .
        (j_mid, αj_mid) = locate(y_mid, g.y)
        u_mid = interpolate(u[i_mid, j_mid], u[i_mid+1, j_mid+1], [αi_mid; αj_mid])

        ## Find the start point s
        x_s, y_s = [g.x[i]; g.y[j]] - Δt*u_mid
        (i_s, αi_s) = locate(x_s, g.x)
        (j_s, αj_s) = locate(y_s, g.y)
        q_[i,j] = interpolate(q[i_s; j_s], q[i_s+1; j_s+1], [αi_s; αj_s])
    end
    return q_
end
