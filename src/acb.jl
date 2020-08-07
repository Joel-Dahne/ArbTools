function real_abs(x::acb; analytic = false)
    y = parent(x)()
    ccall(("acb_real_abs", Nemo.libarb), Cvoid, (Ref{acb}, Ref{acb}, Cint, Clong),
          y, x, analytic, parent(x).prec)
    return y
end

"""
    rel_accuracy_bits(x)
> Compute the relatively accuracy of the ball `x` in bits.
"""
rel_accuracy_bits(x::acb) = ccall(("arb_rel_accuracy_bits", Nemo.libarb), Int,
                                  (Ref{acb},), x)
