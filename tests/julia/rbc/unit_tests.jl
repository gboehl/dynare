include("../../../julia/Utils.jl")
include("rbc1Dynamic.jl")


function test_firstDerivatives()
    endo_nbr = 6
    exo_nbr = 1

    z = [3.473708662732661,2.1382117680737565e-49,3.473708662732661,1.1003308743750966,0.4458199947535853,1.0308567011204435,1.5,2.1382117680737565e-49,1.1003308743750966,0.4458199947535853,1.0308567011204435,1.5]
    params = [0.99,0.357,30.0,0.45,0.5,0.02,0.95,1.5,0.01]
    steady_state = [3.473708662732661,1.1003308743750966,0.4458199947535853,1.0308567011204435,1.5,2.1382117680737565e-49]

    jacobian1 = zeros(Float64, endo_nbr, length(z)+exo_nbr)
    residuals = zeros(Float64, endo_nbr)

    rbc1Dynamic.dynamic!(z, zeros(1,exo_nbr), params, steady_state, 2, residuals, jacobian1)
    @time rbc1Dynamic.dynamic!(z, zeros(1,exo_nbr), params, steady_state, 2, residuals, jacobian1)
    
    g1 = sparse([1],[1],[1.0],6,19)
    g1.colptr = ones(20)
    g1.rowval = ones(24)
    g1.nzval = ones(24)
    rbc1Dynamic.first_derivatives!(z, zeros(1,exo_nbr), params, steady_state, 2, g1)
    @time rbc1Dynamic.first_derivatives!(z, zeros(1,exo_nbr), params, steady_state, 2, g1)

    @assert norm(jacobian1 - g1[:,[1, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 19]]) < 1e-16
end

@time test_firstDerivatives()
