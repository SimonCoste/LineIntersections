module LineIntersections

using DataStructures, Plots

include("line_tools.jl")
include("generator.jl")
include("bentley_ottman.jl")
# Write your package code here.

export Line
export point_in_line, intersection_with_circle, get_endpoints, ∩
export hyperplanes_poisson, hyperplanes_triangle
export intersections_bentley_ottman, intersections_naive
end