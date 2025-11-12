using Plots

function straight_line_path(p, λ=1000)
    """
    Function that generates a path for SMC.jl
        param p₁: starting point
        param p₂: ending point
        param λ: number of points in the path / granularity
        return: array of points in the path
    """
    x₁, y₁ = p[1]
    x₂, y₂ = p[2]
    x = LinRange(x₁, x₂, λ)
    y = LinRange(y₁, y₂, λ)
    path = [collect(x) collect(y)]
    return path
end

function circular_path(p, R, λ=1000)
    """
    Function that generates a circular path for SMC.jl
        param p₁: starting point
        param p₂: ending point
        param R: radius of the circle
        return: array of points in the path
    """
    center_x = p[1]
    center_y = p[2]
    θ = LinRange(0, 2π, λ)
    x = center_x .+ R * cos.(θ)
    y = center_y .+ R * sin.(θ)
    path = [x, y]
    return path
end

function plot_path_dist(p, path, idx_min, n=1000)
    
    diff_path = [p .+ t .* (path[idx_min] .- p) for t in LinRange(0, 1, n)]
    x = [diff_path[i][1] for i in 1:length(diff_path)]
    y = [diff_path[i][2] for i in 1:length(diff_path)]
    return x, y
end

