module LineIntersections

using DataStructures, Plots, Distributions

include("line_tools.jl")
include("artist.jl")
include("generator.jl")
include("bentley_ottman.jl")

export Line
export point_in_line, intersection_with_circle, get_endpoints, ∩
export hyperplanes_poisson, hyperplanes_triangle
export find_all_intersections, intersections_naive
export draw_lines
end