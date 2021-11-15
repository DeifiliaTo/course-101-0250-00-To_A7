using Test, ReferenceTests, BSON

include("../scripts/diffusion_nl_1D.jl")

include("../scripts/viscous_NS_2D.jl")

## Reference Tests with ReferenceTests.jl
# We put both arrays X and T into a BSON.jl and then compare them

"Compare all dict entries"
comp(d1, d2) = keys(d1) == keys(d2) &&
    all([ v1≈v2 for (v1,v2) in zip(values(d1), values(d2))])

# run the model
#H, qx = diffusion_1D()
X, P  = acoustic_2D()

# Test just at some random indices. As for larger models,
# storing the full output array would create really large files!
inds = Int.(ceil.(LinRange(1, length(X), 12)))

d = Dict(:X=> X[inds], :P=>P[inds])
@testset "Ref-tests" begin
#    @test_reference "reftest-files/H.bson" d by=comp
#    @test_reference "reftest-files/qx.bson" d by=comp
     @test_reference "reftest-files/X2.bson" d by=comp
end