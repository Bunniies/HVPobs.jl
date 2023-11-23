module Fits

using LsqFit, LinearAlgebra, ForwardDiff, Statistics, Optim
using ADerrors, LaTeXStrings

struct FitRes
    dof::Int64
    param::Vector{uwreal}
    chi2::Float64
    chi2exp::Float64
    pval::Union{Float64, Nothing}
    
    function FitRes(dof, param, chi2, chi2exp, pval::Union{Float64, Nothing}=nothing)
        
        return new(dof, param, chi2, chi2exp, pval)
    end
end
function Base.show(io::IO, res::FitRes)
    println(io, "FitRes:")
    println(io, " - χ2/χ2exp: " , res.chi2,"/", res.chi2exp, " (", res.chi2/res.chi2exp,")") 
    println(io, " - dof:      ", res.dof)   
    println(io, " - pval:     ", res.pval)
end
export FitRes

include("FitPval.jl")
export get_pvalue

include("FitRoutines.jl")
export fit_routine

include("FitChisq.jl")

end