using Convex
#using ECOS
#using SCS
#using COSMO
#using Tulip
using Random
using LaTeXStrings
using Statistics
using Plots
using LaTeXStrings
using Clarabel
import MathOptInterface as MOI

# This code is a direct implementation of the algorithm described in
# [JV2020] Jiang and Vavasis, ``Certifying clusters from sum-of-norms clustering''

# Build and solve the convex optimization problem
# for several choices of lambdas
function doOptimization(pts,lambdas)
    N = length(pts)
    d = length(pts[1])
    xx = []
    yy = []
    zz = []
    ss = []
    tt = []
    uu = []
    for i = 1:N
        push!(xx,Variable(d))
        push!(zz,Variable(d))
        push!(uu,Variable(1))
        push!(ss,Variable(1))
        for j = (i+1):N
            push!(yy,Variable(d))
            push!(tt,Variable(1))
        end
    end
    lambda = Variable(1)
    #lambda = 1
    sumsquarepart = 1/N*sum(ss)
    closetogetherpart = 1/N^2*lambda*sum(tt)
    constraints = []
    delta_idxs = []
    for i = 1:N
        push!(constraints,1/N*(xx[i]-zz[i]-pts[i]) == 0)
        push!(constraints,1/N*(ss[i]-uu[i]-1) == 0)
        push!(constraints,ss[i] >= norm(vcat(zz[i],uu[i]),2)) # hcat trick from https://discourse.julialang.org/t/convex-jl-norm-of-multiple-variables/98165
    end
    k = 0
    for i = 1:N
        for j = (i+1):N
            k += 1
            push!(constraints,1/N^2*(xx[i]-xx[j]-yy[k]) == 0)
            push!(delta_idxs,length(constraints))
            push!(constraints,tt[k] >= norm(yy[k]))
        end
    end
    problem = minimize(sumsquarepart+closetogetherpart,constraints...)
    optimizer = MOI.OptimizerWithAttributes(Clarabel.Optimizer)
    MOI.set(optimizer,MOI.RawOptimizerAttribute("tol_gap_abs"),1e-15)
    MOI.set(optimizer,MOI.RawOptimizerAttribute("tol_gap_rel"),1e-15)
    rv = []
    for l in lambdas
        @show l
        fix!(lambda,l)
        solve!(problem,optimizer,warmstart=true)
        xs = transpose(hcat([evaluate(x) for x in xx]...))
        deltas = []
        for k in delta_idxs
            push!(deltas,vec(problem.constraints[k].dual))
        end
        dualitygap = MOI.get(problem.model, MOI.ObjectiveValue()) - MOI.get(problem.model, MOI.DualObjectiveValue())
        push!(rv,(xs,deltas))
    end
    return rv
end

# Round the convex optimization problem so that both the primal and dual solutions are feasible
# When we compute the duality gap, we are subtracting two numbers of roughly
# equal magnitude. This is bad news from a numerical stability standpoint.
# Thus, when we do these rounding computations, we do everything using BigFloats
# to get extra precision.
function roundToFeasibilityHighPrec(pts,lambda,xs,deltas)
    pts = [BigFloat.(pt) for pt in pts]
    deltas = [BigFloat.(d) for d in deltas]
    xs = [BigFloat.(x) for x in xs]
    lambda = BigFloat(lambda)
    N = size(xs)[1]
    d = length(xs[1])
    for i = 1:length(deltas)
        if norm(deltas[i]) > lambda
            deltas[i] = lambda/norm(deltas[i]) * deltas[i]
        end
    end
    ys = []
    zs = []
    ss = []
    us = []
    ts = []
    betas = [BigFloat.(zeros(d)) for i = 1:N]
    gammas = []
    k = 0
    for i = 1:N
        for j = (i+1):N
            k += 1
            push!(ys,copy(xs[i]-xs[j]))
            push!(ts,norm(ys[k]))
            betas[i] = betas[i] - 1/N .* deltas[k]
            betas[j] = betas[j] + 1/N .* deltas[k]
        end
        push!(zs,copy(xs[i]-pts[i]))
        push!(ss,1/2*(1+norm(zs[i])^2))
        push!(us,1/2*(-1+norm(zs[i])^2))
    end
    for i = 1:N
        push!(gammas,1/2*(1-norm(betas[i])^2))
    end
    primal_value = 1/N*sum(ss) + lambda*1/N^2*sum(ts) # 
    dual_value = 1/N * sum(gammas) + 1/N * sum(dot.(pts,betas))
    duality_gap = primal_value - dual_value
    #@show duality_gap^(3/4)
    @show duality_gap
    @assert duality_gap >= 0
    return (deltas,primal_value,dual_value,ss,gammas,zs,betas,us,ys)
end

# build the candidate clusters in accordance with the scheme described on [JV2020, p. 8]
function buildCandidateClusters(xs,deltas,duality_gap)
    N = length(xs)
    rv = zeros(Int64,N)
    k = 0
    while true
        k += 1
        j = findfirst(x -> x == 0,rv)
        if j == nothing
            break
        end
        for l = j:N
            if norm(xs[j] - xs[l])<=duality_gap^(3/4)
                rv[l] = k
            end
        end
    end
    return rv
end

# the notation delta_<i,j> defined in [JV2020, (4)]
function bracketik(N,ls,i,j)
    if i > j
        return -bracketik(N,ls,j,i)
    else
        k = 0
        for ii = 1:N
            for jj = (ii+1):N
                k += 1
                if (ii,jj) == (i,j)
                    return ls[k]
                end
            end
        end
    end
    return 0
end

# check the CGR subgradient conditioned defined at the top of [JV2020, p. 9]
function checkCGRCondition(pts,deltas,xs,ss,gammas,zs,betas,us,clusters,lambda)
    N = length(pts)
    sigma1s = ss .* (1 .- gammas) + dot.(zs,betas) + us .* gammas
    sigma2s = ss .* betas + (1 .- gammas) .* zs
    sigma3s = ss .* gammas + (1 .- gammas) .* us
    omegas = sigma3s ./ ss .* zs + 1 ./ ss .* sigma2s
    qs = []
    k = 0
    rprimes = []
    for i = 1:N
        rprime = 1/N * sum(clusters .== clusters[i])
        push!(rprimes,rprime)
        for j = (i+1):N
            k += 1
            if clusters[i] == clusters[j]
                qij = -deltas[k]+1/rprime*(xs[i]-xs[j] - omegas[i] + omegas[j])
                for l = 1:N
                    if clusters[l] != clusters[i]
                        qij .-= 1/(N*rprime) * (bracketik(N,deltas,i,l) - bracketik(N,deltas,j,l))
                    end
                end
                push!(qs,qij)
            else
                push!(qs,0)
            end
        end
    end
    return sum(norm.(qs) .>= lambda) == 0
end

function checkSeparationCondition(clusters,xs,duality_gap)
    N = length(xs)
    nclusters = maximum(clusters)
    for I = 1:nclusters
        for J = (I+1):nclusters
            xbar = mean([xs[i] for i in 1:length(xs) if clusters[i] == I || clusters[i] == J])
            D = 0
            for l = 1:N
                if clusters[l] == I || clusters[l] == J
                    D += 1/N * norm(xs[l] - xbar)^2
                end
            end
            if D <= 2*duality_gap
                @show "Separation condition failed"
                @show D,duality_gap
                return false
            end
        end
    end
    return true
end

function certifyClusters(pts,clusters,xs,deltas,duality_gap,ss,gammas,zs,betas,us,lambda)
    return checkSeparationCondition(clusters,xs,duality_gap) && checkCGRCondition(pts,deltas,xs,ss,gammas,zs,betas,us,clusters,lambda)
end

function doClustering(pts,lambdas)
    optimization_results = doOptimization(pts,lambdas)
    rv = []
    for ((xs, deltas),lambda) in zip(optimization_results,lambdas)
        @show lambda
        xs = [xs[i,:] for i = 1:size(xs,1)]
        deltas,primal_value,dual_value,ss,gammas,zs,betas,us,ys = roundToFeasibilityHighPrec(pts,lambda,xs,deltas)
        duality_gap = primal_value - dual_value
        @show duality_gap
        if duality_gap < 0
            println("duality gap is negative??")
            @assert duality_gap >= 0
        end
        clusters = buildCandidateClusters(xs,deltas,duality_gap)
        if certifyClusters(pts,clusters,xs,deltas,duality_gap,ss,gammas,zs,betas,us,lambda)
            @show clusters
            push!(rv,(clusters,xs))
        else
            println("Certification of clusters failed")
            push!(rv,Nothing)
        end
    end
    return rv
end


function plotClusters(pts,xs,clusters ; show_cluster_reps = false, size=(290,290))
    nclusters = maximum(clusters)
    @show nclusters
    @show clusters
    mycolormap = distinguishable_colors(nclusters, [RGB(1,1,1), RGB(0,0,0)], dropseed=true,transform=deuteranopic) # color-blind friendly, supposedly
    display(mycolormap)
    mycolors = [mycolormap[i] for i in clusters]
    default(
      linewidth=1, 
      framestyle=:box, 
      label=nothing, 
      grid=false
    )
    plt = scatter([pt[1] for pt in pts], [pt[2] for pt in pts],
                  color=mycolors,
                  markerstrokecolor=mycolors,
                  markersize=3,
                  aspect_ratio=:equal,
                  dpi = 100,
                  size=size, #(290,115),
                  legend=false)
    if show_cluster_reps
        scatter!(plt,[x[1] for x in xs], [x[2] for x in xs],
                 aspect_ratio=:equal,
                 color=mycolors,
                 markersize=4,
                 markershape=:+,
                 legend=false)
    end
    return plt
end


# Examples

function stochasticBallModel(d,n)
  rv = []
  for j = 1:n
    a = ones(d)
    while norm(a) >= 1
      a = 2*rand(d)-ones(d)
    end
    push!(rv,a)
  end
  return rv
end

function twoStochasticBalls(d,n,r)
  n1 = sum(rand([0,1],n))
  n2 = n - n1
  ball1 = stochasticBallModel(d,n1)
  ball2 = stochasticBallModel(d,n2)
  for b in ball1
    b[1] += r
  end
  for b in ball2
    b[1] -= r
  end
  return vcat(ball1, ball2), n1
end
    
    
function stochasticSphereModel(d,n)
  rv = []
  for j = 1:n
    a = 2*rand(d)-ones(d)
    a = a ./ norm(a)
    push!(rv,a)
  end
  return rv
end

function twoStochasticSpheres(d,n,r)
  n1 = sum(rand([0,1],n))
  n2 = n - n1
  ball1 = stochasticSphereModel(d,n1)
  ball2 = stochasticSphereModel(d,n2)
  for b in ball1
    b[1] += r
  end
  for b in ball2
    b[1] -= r
  end
  return vcat(ball1, ball2), n1
end


function evenlySpacedCircle(n)
  rv = []
  for j = 1:n
    push!(rv,[cos(2*pi*j/n),sin(2*pi*j/n)])
  end
  return rv
end

function twoEvenlySpacedCircles(n,r)
  n1 = Int(floor(n/2))
  n2 = Int(ceil(n/2))
  ball1 = evenlySpacedCircle(n1)
  ball2 = evenlySpacedCircle(n2)
  for b in ball1
    b[1] += r
  end
  for b in ball2
    b[1] -= r
  end
  return shuffle(vcat(ball1, ball2)), n1
end


function threeCircles(n,r1,r2,r3,s1,s2,s3)
  n1 = Int(floor(n/3))
  n2 = Int(ceil(n/3))
  n3 = n-n1-n2
  ball1 = evenlySpacedCircle(n1)
  ball2 = evenlySpacedCircle(n2)
  ball3 = evenlySpacedCircle(n2)
  for b in ball1
      b[1] *= s1
      b[2] *= s1
  end
  for b in ball2
      b[1] *= s2
      b[2] *= s2
  end
  for b in ball3
      b[1] *= s3
      b[2] *= s3
  end
  for b in ball1
    b[1] += r1
  end
  for b in ball2
    b[1] += r2
    b[2] += r3
  end
  return vcat(ball1, ball2,ball3)
end

function evenlySpacedBall2d(delta)
    pts = []
    for i=-1:delta:1
        for j = -1:delta:1
            if norm([i,j])<=1
                push!(pts,[i,j])
            end
        end
    end
    return pts
end

function evenlySpacedRoundBall2d(delta)
    pts = [[0.0,0.0]]
    for r = delta:delta:1
        for theta = 0:pi/(floor(pi/(delta/r))):2pi
            push!(pts,[r*cos(theta),r*sin(theta)])
        end
    end
    return pts
end


function sizes(l)
    N = length(l)
    sizes = []
    for i = 1:N
        sz = 0
        for j = 1:N
            if l[j] == l[i]
                sz += 1
            end
        end
        push!(sizes,sz)
    end
    return sizes
end
