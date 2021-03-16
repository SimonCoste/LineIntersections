function hyperplanes_poisson(R, n::Int)::Array{Line}
    """
    Inputs:
        max radius R 
        number of lines n

    Outputs:
        n lines with parameters (θ, r) uniformly random 
        in the set [0, 2π)×[0,R).
    """
    α = (2*pi) .* rand(n) ; radii = R * rand(n)
    return Array([Line(α[k], radii[k]) for k in 1:n])
end

function hyperplanes_regrad(R, n::Int)::Array{Line}
    """
    Inputs:
        max radius R 
        number of lines n

    Outputs:
        n lines with parameters (θ, r) uniformly random 
        in the set [0, 2π)×[0,R).
    """
    α = (2*pi) .* rand(n) ; radii = Array([R*(i/n) for i in 1:n])
    return Array([Line(α[k], radii[k]) for k in 1:n])
end

function hyperplanes_triangle(R, n::Int)::Array{Line}
    """
    Inputs:
        max radius R 
        number of lines n

    Outputs:

    """
    radii = R * rand(n)
    angles = [pi/2, pi/2 + 2*pi/3, pi/2 + 4*pi/3]
    α = rand(angles ∪ -angles, n) 
    return Array([Line(α[k], radii[k]) for k in 1:n])
end

#hyperplanes_triangle(1, 10)

