@testset "line tools" begin
    
    line = Line(pi/2, 1)
    @test point_in_line(line, 0.)==1.
end