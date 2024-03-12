function Base.abs(uw::uwreal)
    sgn = sign(value(uw))
    return sgn * uw 
end