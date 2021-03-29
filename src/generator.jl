using Distributions

function hyperplanes_poisson(R, λ)::Array{Line}
    """
    Inputs:
        max radius R 
        intensity λ

    Outputs:
        N lines with parameters (θ, r) uniformly random, 
        in the set [0, 2π)×[0,R), where N=Poi(2πRλ) .
    """
    n = rand(Poisson(2*pi*R*λ), 1)[1]
    α = (2*pi) .* rand(n) ; radii = R * rand(n)
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

function hyperplanes_anisotropic_poisson(R, K, Λ)::Array{Line}
    """
    Inputs:
        - radius R 
        - 1d array of directions K
        - 1d array of intensities Lambda with same dimension as k
    """
    @assert length(K)==length(Λ) throw(DimensionMismatch)

    n = rand(Poisson(2*pi*R*Λ[1]), 1)[1]
    output = [Line(K[1], x) for x in 2*R*rand(n).-R]
    if length(K)==1
        return output
    else
        for i in 2:length(K)
            n = rand(Poisson(2*pi*R*Λ[i]), 1)[1]
            newdir = [Line(K[i], x) for x in 2*R*rand(n).-R]
            output = output ∪ newdir 
        end
        return output
    end
end