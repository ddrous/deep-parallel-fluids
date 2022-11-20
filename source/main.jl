## Imports
include("solver.jl")
include("visualization.jl")

# exit()

## Initiate a NavierStokes problem
function initnavierstokes(g::Geometry; ρ=1.0, fval::Float64=-1.0)
    u = zeros(Float64, (g.Nx+1, g.Ny+1, 2))   ## Anticipating zero velocitity in solid
    p = ones(Float64, (g.Nx+1, g.Ny+1))


    u[3:g.Nx-1, 3:g.Ny-1, :] .= 1.0     ## Non-zero velocity in legining of fluid

    p[:, g.Ny:g.Ny+1] .= 0.0   ## Pressure is zero in air at east

    f = fill(fval, (g.Nx+1, g.Ny+1, 2))
    f[:, :, 1] .= 0.0

    return NavierStokes(u, p, f, ρ)
end



### Run the code in a function: no global variables !
function main()
    geo = Geometry(; Lx=1, Ly=1, Nx=50, Ny=30)
    pb = initnavierstokes(geo; ρ=1.0, fval=-1.0)

    Δt = 5 * min(geo.Δx, geo.Δy) / maximum(pb.u)
    @time solvenavierstokes(pb, Δt, geo)
    # visualizenavierstokes(pb, geo)
    makenavierstokesgif(pb, geo)
end


main()

# exit()
