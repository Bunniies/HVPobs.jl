module LinAlg

using ADerrors
import ADerrors.uwerr
# import Base.:*

# struct for Rational type with uwreal objects
mutable struct Rationaluw <: Real
    num::uwreal
    den::uwreal
    function Rationaluw(num::uwreal, den::uwreal)
        norm = gcd(Int(value(num)), Int(value(den))) 
        println(norm)
        num /= norm
        den /= norm
        return new(num, den)
    end
end
function Base.show(io::IO, aa::Rationaluw)
    println(io, aa.num,"//",aa.den)
end
function ADerrors.uwerr(aa::Rationaluw)
    uwerr(aa.num)
    uwerr(aa.den)
end
export Rationaluw, uwerr

include("LinAlgBroadcast.jl")

end