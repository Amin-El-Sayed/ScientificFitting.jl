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

"""Fixed-point rendering of `x` with `digits` decimal places (0 for integers)."""
_fixed_decimals(x::Real, digits::Int) =
    Printf.format(Printf.Format("%.$(max(digits, 0))f"), round(Float64(x); digits=digits))

"""
    _value_error_strings(value, err; sigdigits=5) -> (String, String)

Round a value/uncertainty pair for reporting: the uncertainty keeps two
significant digits and the value is rounded to the same decimal place, so
`15.5551 +/- 0.10753` renders as `("15.56", "0.11")`. Pairs outside the
fixed-point range (non-finite or non-positive uncertainties, uncertainties
below 1e-6 or at 1e5 and above) fall back to the plain significant-digit
format of `_fmt_value` for both entries.
"""
function _value_error_strings(value::Real, err::Real; sigdigits::Int=5)
    digits = _error_decimals(value, err)
    digits === nothing && return (
        _fmt_value(value; sigdigits=sigdigits), _fmt_value(err; sigdigits=sigdigits))
    return (_fixed_decimals(value, digits), _fixed_decimals(err, digits))
end

"""
    _value_error_strings(value, minus, plus; sigdigits=5) -> (String, String, String)

Asymmetric-uncertainty variant: the smaller uncertainty fixes the decimal
place (two significant digits), and all three numbers render on it.
"""
function _value_error_strings(value::Real, minus::Real, plus::Real; sigdigits::Int=5)
    digits = _error_decimals(value, min(minus, plus))
    fallback_digits = _error_decimals(value, max(minus, plus))
    digits === nothing && (digits = fallback_digits)
    digits === nothing && return (
        _fmt_value(value; sigdigits=sigdigits),
        _fmt_value(minus; sigdigits=sigdigits),
        _fmt_value(plus; sigdigits=sigdigits))
    return (_fixed_decimals(value, digits),
            _fixed_decimals(minus, digits),
            _fixed_decimals(plus, digits))
end

# Decimal place carrying two significant digits of `err`, or `nothing` when
# the pair has no compact fixed-point rendering.
function _error_decimals(value::Real, err::Real)
    (isfinite(value) && isfinite(err) && err > 0) || return nothing
    digits = 1 - floor(Int, log10(abs(err)))
    (-4 <= digits <= 6) || return nothing
    return digits
end
