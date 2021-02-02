function integer_roots(p::Nemo.fmpq_poly)
    result = []
    for (f, _) in factor(p)
        if (degree(f) == 1) && (coeff(f, 1) == 1)
            push!(result, -evaluate(f, 0))
        end
    end
    return result
end
