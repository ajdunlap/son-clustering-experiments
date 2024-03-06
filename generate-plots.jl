include("sonclustering.jl")
using Printf

pythonplot()

begin
    lambdas = 1.09:0.0025:1.125
    meshes = 0.0625:0.0075:0.115

    for mesh in meshes
        pts = evenlySpacedBall2d(mesh)
        results = doClustering(pts,lambdas)
        vals = [(l,maximum(r[1])) for (l,r) in zip(lambdas,results) if r != Nothing]
        ls = [v[1] for v in vals]
        rs = [v[2] for v in vals]
        plt = scatter(ls,rs,xlabel=L"$\lambda$",ylabel="no. of clusters",legend=false,size=(290,140),xlim=(1.09,1.125),xticks=[1.09,1.10,1.11,1.12],yticks=[1,length(pts)])
        savefig(plt,@sprintf("find-lambda1-%.4f.tex",mesh))
    end

    pts = threeCircles(100,3.5,1.5,3,0.7,0.5,0.6)

    clusters,xs = doClustering(pts,1.1)[1]
    plt = plotClusters(pts,xs,clusters;show_cluster_reps=true)
    savefig(plt,"three-circles-1.4.tex")

    clusters,xs = doClustering(pts,2.4)[1]
    plt = plotClusters(pts,xs,clusters;show_cluster_reps=true)
    savefig(plt,"three-circles-2.4.tex")

    clusters,xs = doClustering(pts,3.0)[1]
    plt = plotClusters(pts,xs,clusters;show_cluster_reps=true)
    savefig(plt,"three-circles-3.0.tex")

    clusters,xs = doClustering(pts,3.4)[1]
    plt = plotClusters(pts,xs,clusters;show_cluster_reps=true)
    savefig(plt,"three-circles-3.4.tex")

    clusters,xs = doClustering(pts,3.6)[1]
    plt = plotClusters(pts,xs,clusters;show_cluster_reps=true)
    savefig(plt,"three-circles-3.6.tex")

    include("200stochasticballs.jl")
    clusters,xs = doClustering(pts,2)[1]
    plt = plotClusters(pts,xs,clusters;show_cluster_reps=false,size=(290,120))
    plot!(plt,t->cos(t)-1.05,sin,0,2pi,linestyle=:dot,color=:grey)
    plot!(plt,t->cos(t)+1.05,sin,0,2pi,linestyle=:dot,color=:grey)
    savefig(plt,"shattered-balls.tex")
    clusters,xs = doClustering(pts,2.15)[1]
    plt = plotClusters(pts,xs,clusters;show_cluster_reps=false,size=(290,120))
    plot!(plt,t->cos(t)-1.05,sin,0,2pi,linestyle=:dot,color=:grey)
    plot!(plt,t->cos(t)+1.05,sin,0,2pi,linestyle=:dot,color=:grey)
    savefig(plt,"cohesive-balls.tex")
end

begin
    N = 8

    pts,n1 = twoEvenlySpacedCircles(2*N,1.7)

    lambdas = 2 .* [1.55,1.65,1.75]
    results = doClustering(pts,lambdas)
    for ((clusters,xs),l) in zip(results,lambdas)
        plt = plotClusters(pts,xs,clusters;show_cluster_reps=false,size=(290,120))
        savefig(plt,@sprintf("two-octagons-%.2f.tex",l))
    end
end
