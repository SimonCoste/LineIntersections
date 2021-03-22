@testset "algorithms" begin
    
    number_of_lines = 10
    H1 = hyperplanes_poisson(1, number_of_lines)
    res1a = find_all_intersections(H1)
    res1b = intersections_naive(H1)

    @test length(res1a)==length(res1b)

end