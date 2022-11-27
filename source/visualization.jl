
include("structures.jl")
include("geometry.jl")
include("solver.jl")

using Plots, SparseArrays, LinearAlgebra
# plotlyjs()


## visualize a time-step of the navierstokes problem's solution
function visualizenavierstokes(l::NavierStokes, g::Geometry)
    ## Display solution in a 2D plane
    println("\n----Solution----")
    display(matto2d(l.u))
    println()

    ## Construct plot data: y axis first, then x axis (Julia convention)
    x=range(0,g.Lx,length=g.Nx+1)
    y=range(0,g.Ly,length=g.Ny+1)
    data = [norm(l.u[i,j,:]) for j in 1:g.Ny+1, i in 1:g.Nx+1]

    anns = [(0, 5, Plots.text("  N=$(g.Nx)")),
            (0, 4.7, Plots.text("M=$(g.Ny)")),
            (1, 5, Plots.text("Δx=$(g.Lx/g.Nx)")),
            (1, 4.7, Plots.text("Δy=$(g.Ly/g.Ny)")),
            (2, 5, Plots.text("f=$(l.f[1,1])")),
            (2, 4.7, Plots.text("   bc=$(l.bc)")),
            (0, 0, Plots.text("RDN")),
            (2, 5.6, Plots.text(" "))
            ]
    title = "Numerical solution of the Navier-Stokes problem"

    plot(x, y, data,
        xlabel="x", ylabel="y", zlabel="u", 
        camera=(30, 30), 
        st=:surface,
        annotations=anns,
        title=title,
        titlefontsize=18)
end


## Visiualize a gif of the NavierStokes problem's solution
function makenavierstokesgif(l::NavierStokes, g::Geometry)
    x=range(0, g.Lx, length=g.Nx+1)
    y=range(0, g.Ly, length=g.Ny+1)

    anns = [(0, 5, Plots.text("  N=$(g.Nx)")),
            (0, 4.7, Plots.text("M=$(g.Ny)")),
            (1, 5, Plots.text("Δx=$(g.Lx/g.Nx)")),
            (1, 4.7, Plots.text("Δy=$(g.Ly/g.Ny)")),
            (2, 5, Plots.text("f=$(round(l.f[1,1,2]))")),
            (2, 4.7, Plots.text("   u_wall=$(round(l.u[1,1]))")),
            (0, 0, Plots.text("RDN")),
            (2, 5.6, Plots.text(" "))
            ]
    title = "Numerical solution of the NavierStokes problem"


    # @gif for t in 0:Δt:l.T 
    @gif for t in 1:1:1
        Δt = 0.005 * min(g.Δx, g.Δy) / maximum(l.u)
        solvenavierstokes(l, Δt, g)

        data = [norm(gatherfield_xy(i,j,l.u,l.v)) for j in 1:g.Ny, i in 1:g.Nx]

        # plot(x, y, data,
        # xlabel="x", ylabel="y", zlabel="u", 
        # camera=(30, 30), 
        # st=:surface,
        # annotations=anns,
        # title=title,
        # titlefontsize=18)

        heatmap(x, y, data,
                xlabel="x", ylabel="y", zlabel="u",
                # annotations=anns,
                title=title,
                titlefontsize=12)
    end

end
