using LineIntersections
using Test



@testset "LineIntersections.jl" begin
    include("test_line_tools.jl")
    include("test_bentley_ottman.jl")
end
