function real_abs(x::acb; analytic = false)
    # FIXME: The version in Arb contains a bug in one case, we check for this
    # manually. Once this is fixed in Arb we can remove it from here.
    if !analytic && !isnonnegative(real(x)) && !isnegative(real(x))
        return return parent(x)(setunion(real(x), -real(x)), setunion(imag(x), -imag(x)))
    end
    y = parent(x)()
    ccall(("acb_real_abs", Nemo.libarb), Cvoid, (Ref{acb}, Ref{acb}, Cint, Clong),
          y, x, analytic, precision(parent(x)))
    return y
end

"""
    rel_accuracy_bits(x)
> Compute the relatively accuracy of the ball `x` in bits.
"""
rel_accuracy_bits(x::acb) = ccall(("arb_rel_accuracy_bits", Nemo.libarb), Int,
                                  (Ref{acb},), x)
