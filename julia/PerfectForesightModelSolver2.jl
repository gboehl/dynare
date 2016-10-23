module PerfectForesightModelSolver2

##
 # Copyright (C) 2016 Dynare Team
 #
 # This file is part of Dynare.
 #
 # Dynare is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # Dynare is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
##

import DynareModel.Model
import DynareOutput.Output
import DynareOptions.Options

export simulate_perfect_foresight_model!

function simulate_perfect_foresight_model!(endogenousvariables::Matrix{Float64}, exogenousvariables::Matrix{Float64}, steadystate::Vector{Float64}, model::Model, options::Options)

    lead_lag_incidence = model.lead_lag_incidence

    nyp = countnz(lead_lag_incidence[1,:])
    ny0 = countnz(lead_lag_incidence[2,:])
    nyf = countnz(lead_lag_incidence[3,:])

    ny = length(model.endo)
    nd = nyp+ny0+nyf

    periods = options.pfmsolver.periods
    params = model.params

    tmp = lead_lag_incidence[2:3,:]'
    i_cols_A1 = find(tmp)
    i_cols_1  = tmp[i_cols_A1]

    tmp = lead_lag_incidence[1:2,:]'
    i_cols_AT = find(tmp)
    i_cols_T  = tmp[i_cols_AT]

    tmp = lead_lag_incidence[2,:]'
    i_cols_A0 = find(tmp)
    i_cols_0  = tmp[i_cols_A0]

    i_cols_j = collect(1:nd)
    i_upd = ny+collect(1:periods*ny)

    Y = vec(endogenousvariables)
    z = Y[find(lead_lag_incidence')]

    jacobian = zeros(Float64, ny, length(z)+length(model.exo))
    residuals = zeros(Float64, ny)

    println("\nMODEL SIMULATION:\n")

    rd = zeros(Float64, periods*ny)
    A = spzeros(periods*ny,periods*ny)
    convergence = false
    iteration = 0

    ws = PerfectForesightModelWS(periods,model)
    
    while !convergence
        iteration += 1
        perfect_foresight_model!(Y,exogenousvariables,model,steadystate,ws,rd)
        err = maximum(abs(rd))
        println("Iter. ", iteration, "\t err. ", round(err, 12))
        if err<options.pfmsolver.tolf
            iteration -= 1
            convergence = true
            break
        end
        A = perfect_foresight_model!(Y,exogenousvariables,model,steadystate,ws,A)
        @time dy = A\rd
        Y[i_upd] -= dy
        if maximum(abs(dy))<options.pfmsolver.tolx
            convergence = true
        end
    end
    if convergence
        println("\nPFM solver converged in ", iteration, " iterations!\n")
        endogenousvariables = reshape(Y, ny, periods+2)
    end
end

type PerfectForesightModelWS
    periods::Int64
    iA::Array{Int64,1} 
    jA::Array{Int64,1} 
    vA::Array{Float64,1} 
    iP::Array{Int64,1} 
    jP::Array{Int64,1} 
    vP::Array{Float64,1} 
    function PerfectForesightModelWS(n,model)
        periods = n
        nzd = model.nnzderivatives[1]
        iA = zeros(Int64,n*nzd)
        jA = zeros(Int64,n*nzd)
        vA = zeros(Float64,n*nzd)
        iP = zeros(Int64,nzd)
        jP = zeros(Int64,nzd)
        vP = zeros(Float64,nzd)
        new(periods,iA,jA,vA,iP,jP,vP)
    end
end

function perfect_foresight_model!(Y::Array{Float64,1},exogenousvariables::Array{Float64,2},model::Model,steadystate::Array{Float64,1},ws::PerfectForesightModelWS,residuals::Array{Float64,1})
    periods = ws.periods
    ny = length(model.endo)
    i_rows = collect(1:ny)
    i_cols = find(model.lead_lag_incidence')
    m = 0
    for it = 2:(periods+1)
        res = sub(residuals,m+1:m+model.eq_nbr)
        Yview = sub(Y,i_cols)
        model.dynamic(Yview, exogenousvariables, model.params, steadystate, it, res)
        m += model.eq_nbr
        i_cols += ny
    end    
end

function perfect_foresight_model!(Y::Array{Float64,1},exogenousvariables::Array{Float64,2},model::Model,steadystate::Array{Float64,1},ws::PerfectForesightModelWS,A::SparseMatrixCSC{Float64,Int64})
    periods = ws.periods
    ny = length(model.endo)
    i_rows = collect(1:ny)
    i_cols = find(model.lead_lag_incidence')
    m = 0
    offset_r = 0
    offset_c = -ny
    nzd = 0
    for it = 2:(periods+1)
        Yview = sub(Y,i_cols)
        if it == 2
            model.first_derivatives(Y[i_cols], exogenousvariables, model.params, steadystate, it, ws.iP, ws.jP, ws.vP)
            nzd = count(i->(0 .< i .<= 3*ny),ws.jP)
            k1 = 1
            for k = 1:nzd
                if ws.jP[k] > ny
                    ws.iA[k1] = ws.iP[k]
                    ws.jA[k1] = ws.jP[k] + offset_c
                    ws.vA[k1] = ws.vP[k]
                    k1 += 1
                end
            end
            m = k1 - 1
        elseif it==(periods+1)
            model.first_derivatives(Y[i_cols], exogenousvariables, model.params, steadystate, it, ws.iP, ws.jP, ws.vP)
            for k=1:nzd
                if ws.jP[k] <= 2*ny
                    m += 1
                    ws.iA[m] = ws.iP[k] + offset_r
                    ws.jA[m] = ws.jP[k] + offset_c
                    ws.vA[m] = ws.vP[k]
                end
            end
        else
            I = sub(ws.iA,m+1:m+nzd)
            J = sub(ws.jA,m+1:m+nzd)
            V = sub(ws.vA,m+1:m+nzd)
            model.first_derivatives(Y[i_cols], exogenousvariables, model.params, steadystate, it, I,J,V)
            for k=m+1:m+nzd
                ws.iA[k] += offset_r
                ws.jA[k] += offset_c
            end
            m += nzd        
        end
        offset_r += ny
        offset_c += ny
        i_cols += ny
    end
    A = sparse(ws.iA[1:m], ws.jA[1:m], ws.vA[1:m])
    assert(size(A)==(ny*periods,ny*periods))
    return A
end

end
