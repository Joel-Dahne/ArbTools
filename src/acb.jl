function real_abs(x::acb; analytic = false)
    y = parent(x)()
    ccall(("acb_real_abs", Nemo.libarb), Cvoid, (Ref{acb}, Ref{acb}, Cint, Clong),
          y, x, analytic, parent(x).prec)
    return y
end
