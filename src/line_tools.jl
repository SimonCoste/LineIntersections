
import Base: ∩


struct Line
    """
    Lines are described by the polar coordinates (θ, r) 
    of the orthogonal projection of the origin. 
        
    - θ is in [0, 2pi).
        - r is a nonnegative number. 
    """
    θ::Float64
    r::Float64
end




function point_in_line(ℓ::Line, x::Float64)::Float64
    """
    Outputs the unique y such that (x, y) ∈ ℓ.
    """
    r = ℓ.r ; θ = ℓ.θ
    return r/sin(θ) - x/tan(θ)
end

function point_in_line(l::Line, x::Array{Float64})::Array{Float64}
    """
    Computes L(x) for any array x
    """
    r = l.r ; θ = l.θ
    return r/sin(θ) .- ( x ./ tan(θ) )
end


function intersection_with_circle(l::Line, R::Float64=1.)
    """
    Computes the two intersection of line l with the 
    circle of radius R centered at zero. 

    If there are no such intersections, raises an AssertionError. 
    """
    @assert R>l.r "circle radius R must be greater than r, the line distance to 0"
    α = acos(l.r / R)
    x₁ = Array([R * cos(l.θ + α), R * sin(l.θ + α)])
    x₂ = Array([R * cos(l.θ - α), R * sin(l.θ - α)])
    return Array([x₁, x₂])

end


function get_endpoints(L::Array{Line, 1})
    """
    Given an array of n lines, computes the intersections of each
    line with a given circle of radius 1, then outputs the 2n intersection points. 

    Input:
        an array of n Line
    
    Output:
        V, a 2×2n array of Floats. The vector V[:, 2i-1] is the leftpoint of
        line H[i] and V[:, 2i] is the right endpoint.
    """
    n = length(L) #number of lines
    output = zeros(2, 2 * n) #first line is x-coordinate of point
    
    for i in 1:n

        (a,b) = intersection_with_circle(L[i])
        
        if a[1]<b[1]
            output[:, 2*i-1] = a
            output[:, 2*i] = b
        else
            output[:, 2*i-1] = b
            output[:, 2*i] = a
        end
    end
    return output
end


function ∩(l1::Line, l2::Line)::Array{Float64}
    """
    Input: 
        two nonparallel line objects l1 and l2.

    Output: 
        the cartesian coordinates of their intersection point. 
    """
    #@assert l1.θ!=l2.θ "Lines must not be parallel. "
    u = (l2.r / sin(l2.θ)) - (l1.r / sin(l1.θ))
    d = (1/tan(l2.θ) - 1/tan(l1.θ)) 
    x = u/d
    return Array([x, point_in_line(l1, x)])
end

normof(p) = p[1]^2 + p[2]^2


function draw(l::Line)
    x = collect(range(-1, stop=1, length=100))
    plot!(x, point_in_line(l, x))
end



function draw_lines(H)
    x = collect(range(-1, stop=1, length=3))
    v = collect(range(0, stop=2*pi, length=100))
    dessin = plot(sin.(v), cos.(v), color=:gray)
    for l in H
        plot!(dessin, x, point_in_line(l, x), color=:thistle)
    end
    plot!(dessin, xlims=(-1, 1), ylims=(-1, 1))
end

#dessin = draw_lines(H)



