# example of using the same umfpack symbolic factorization to solve two
# linear systems with the same sparsity pattern. Usefull for nonlinear solvers
# timing and memory use is not representative of a perfect foresight problem
# because the test matrices are not block diagonal

using Base.SparseMatrix.UMFPACK

n=100000
a1 = sprand(n,n,0.00001)+speye(n)
k = find(a1)
a2 = deepcopy(a1)
a2[k] = randn(length(k))
b = randn(n,1)

@time u1 = lufact(a1)
@time res = a1\b
println(res[1:10])
@time res = a2\b
println(res[1:10])
u2 = deepcopy(u1)
@time UMFPACK.umfpack_symbolic!(u1)
@time res = u1\b
println(res[1:10])
u2.symbolic = deepcopy(u1.symbolic)
u2.nzval = deepcopy(a2.nzval)
u2.numeric = C_NULL
@time UMFPACK.umfpack_numeric!(u2)
@time UMFPACK.umfpack_numeric!(u2)
@time res = u2\b
println(res[1:10])
