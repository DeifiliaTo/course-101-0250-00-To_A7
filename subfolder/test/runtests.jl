using Test, ReferenceTests, BSON

include("../scripts/diffusion_nl_1D.jl")

## Unit tests
@testset "av" begin
    @test av([1, 2, 3]) ≈ [1.5, 2.5]
    @test av([-1, -2, -3]) ≈ [-1.5, -2.5]
    @test av([0.001, 0.002, 0.003, 0.006]) ≈ [0.0015, 0.0025, 0.0045]
end


## Reference Tests with ReferenceTests.jl
# We put both arrays X and T into a BSON.jl and then compare them

"Compare all dict entries"
comp(d1, d2) = keys(d1) == keys(d2) &&
    all([ v1≈v2 for (v1,v2) in zip(values(d1), values(d2))])

# run the model
H, qx = diffusion_1D()

# Test just at some random indices. As for larger models,
# storing the full output array would create really large files!
inds = [18, 27, 45, 68, 71, 71, 102, 110, 123]

d = Dict(:H=> H[inds], :qx=>qx[inds])
@testset "Ref-tests" begin
    @test_reference "reftest-files/H.bson" d by=comp
    @test_reference "reftest-files/qx.bson" d by=comp
end