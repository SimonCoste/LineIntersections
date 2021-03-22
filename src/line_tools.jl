
import Base:∩,<

##################################################
# Types for lines 
##################################################

abstract type AbstractLine end

struct Line <: AbstractLine
    """
    Lines are described by the polar coordinates (θ, r) 
    of the orthogonal projection of the origin. 
        
    - θ is in [0, 2pi).
        - r is a nonnegative number. 
    """
    θ::Float64
    r::Float64
end

mutable struct SweepLine <: AbstractLine
    """
    A vertical line, for the BO algorithm.
    """
    θ::Float64
    r::Float64
    SweepLine(r) = new(0, r)
end

struct RootedLine <: AbstractLine
    """
    It is a line, rooted at a certain x-position given by
    the sweep-line. 
    """
    θ::Float64
    r::Float64
    ∅::SweepLine
    RootedLine(l::Line, sl::SweepLine) = new(l.θ, l.r, sl)
end

################################################
# Methods 
################################################


function point_in_line(ℓ::AbstractLine, x::Float64)::Float64
    """
    Outputs the unique y such that (x, y) ∈ ℓ.
    """
    r = ℓ.r ; θ = ℓ.θ
    @assert θ!=0 throw(error("No vertical lines."))
    return r/sin(θ) - x/tan(θ)
end

function <(l1::RootedLine, l2::RootedLine)
    """
    A rooted line is smaller than another if the intersection with
    its sweepline is smaller than the other.
    """ 
    a = point_in_line(l1, l1.∅.r)
    b = point_in_line(l2, l2.∅.r)
    return a<b
end

function point_in_line(l::Line, x::Array{Float64})::Array{Float64}
    """
    Computes L(x) for any array x
    """
    r = l.r ; θ = l.θ
    return r/sin(θ) .- ( x ./ tan(θ) )
end

function intersection_with_circle(l::AbstractLine, R::Real)
    """
    Computes the two intersection of line l with the 
    circle of radius R centered at zero. 

    If there are no such intersections, returns nothing.
    Outputs the two points from left to right. 
    """
    if R <= l.r
        return nothing
    else
        α = acos(l.r / R)
        x₁ = Array([R * cos(l.θ + α), R * sin(l.θ + α)])
        x₂ = Array([R * cos(l.θ - α), R * sin(l.θ - α)])
        return isless(x₁, x₂) ? (x₁, x₂) : (x₂, x₁)
    end

end

function get_endpoints(L::Array{Line, 1}, R::Real)
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
        output[:, 2*i-1] = a
        output[:, 2*i] = b
    end
    return output
end

normof(p) = p[1]^2 + p[2]^2

function ∩(l1::AbstractLine, l2::AbstractLine)::Array{Float64}
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