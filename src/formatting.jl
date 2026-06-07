function _fmt_value(x::Real; sigdigits::Int=5)
    if isnan(x)
        return "NaN"
    elseif isinf(x)
        return signbit(x) ? "-Inf" : "Inf"
    end
    return string(round(Float64(x); sigdigits=sigdigits))
end

function _strip_math_delims(s::AbstractString)
    if startswith(s, "\$") && endswith(s, "\$") && ncodeunits(s) >= 2
        return s[2:(end - 1)]
    end
    return s
end
