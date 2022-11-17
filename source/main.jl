## Imports
include("solver.jl")
include("visualization.jl")



## Initiate a NavierStokes problem
function initnavierstokes(g::Geometry; bc::Float64=1.0, fval::Float64=1.0)
    u = ones(Float64, (g.Nx+1,g.Ny+1))
    setbc(u, g)
    f = fill(fval, (g.Nx+1,g.Ny+1))
    return NavierStokes(bc, u, f)
end



### Run the code in a function: no global variables !
function main()
    geo = Geometry(1,1,100,50)
    pb = initnavierstokes(geo; bc=1.0, fval=2.0)
    @time solvenavierstokes(pb, geo)
    # visualizenavierstokes(pb, geo)
    makenavierstokesgif(pb, geo)
end


main()

