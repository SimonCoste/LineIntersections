@testset "algorithms" begin
    
    number_of_lines = 10
    H1 = hyperplanes_poisson(1, number_of_lines)
    res1a = intersections_bentley_ottman(H1)
    res1b = intersections_naive(H1)
    
    H2 = hyperplanes_triangle(1., number_of_lines)
    res2a = intersections_bentley_ottman(H2)
    res2b = intersections_naive(H2)
    @test isapprox(res1a, res1b, atol=1e-9)
    @test isapprox(res2a, res2b, atol=1e-9)
end