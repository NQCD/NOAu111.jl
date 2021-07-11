using NOAu111
using Test
using JuLIP

@testset "AuAu" begin
    at = bulk(:Au)*3
    @test JuLIP.Testing.fdtest(AuAu(), at)
end
