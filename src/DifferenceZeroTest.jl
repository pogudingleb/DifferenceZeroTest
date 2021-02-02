using Base
using AbstractAlgebra, Nemo

include("integer_roots.jl")
include("utils.jl")

abstract type AbstractSigmaPSDomain <: AbstractAlgebra.Ring end
abstract type ComputablePS{T} <: AbstractAlgebra.RingElem end

MAX_ORD = 15

#------------------------------------------------------------------------------
# Basics for computable power series
#------------------------------------------------------------------------------

function ps_compose(f::AbstractAlgebra.AbsSeriesElem, g::AbstractAlgebra.AbsSeriesElem)
    result = zero(parent(f))
    power = one(parent(f))
    for i in 0:precision(f)
        result += power * coeff(f, i)
        power = power * g
    end
    return result
end

#------------------------------------------------------------------------------

function get_truncation(ps::ComputablePS, R::Generic.SeriesRing)
    prec = max_precision(R)
    if prec > length(cached(ps))
        truncation = get_truncation_raw(ps, R)
        update_cached!(ps, truncation)
        return truncation
    end
    t = gen(R)
    return sum([ps.cached[i] * t^(i - 1) for i in 1:prec])
end

function get_truncation(ps::ComputablePS{T}, n::Integer) where T <: RingElement
    psring, t = PowerSeriesRing(base_ring(ps), n, "t"; model=:capped_absolute)
    return get_truncation(ps, psring)
end

#------------------------------------------------------------------------------

function cached(ps::ComputablePS)
    return ps.cached
end

#------------------------------------------------------------------------------

function update_cached!(ps::ComputablePS, trunc::AbsSeriesElem)
    ps.cached = [coeff(trunc, i) for i in 0:(precision(trunc) - 1)]
end

#------------------------------------------------------------------------------

function AbstractAlgebra.base_ring(ps::ComputablePS)
    return ps.base_ring
end

#------------------------------------------------------------------------------

function valuation(fs::Array{<: ComputablePS{T}, 1}) where T <: RingElement
    prec = 5
    while true
        R, t = PowerSeriesRing(base_ring(fs[1]), prec, "t"; model=:capped_absolute)
        vals = map(f -> AbstractAlgebra.valuation(get_truncation(f, R)), fs)
        if any(map(v -> v != prec, vals))
            return (min(vals...), findfirst(v -> v == min(vals...), vals))
        end
        prec = 2 * prec
    end
end

function valuation(f::ComputablePS)
    return valuation([f])[1]
end

#------------------------------------------------------------------------------
# Particular computable power series used in the examples
#------------------------------------------------------------------------------

struct ConstPS{T} <: ComputablePS{T}
    base_ring::AbstractAlgebra.Ring
    c::T

    function ConstPS{T}(a::RingElement) where T <: RingElement
        br = parent_type(T)()
        return new(br, br(a))
    end

    function ConstPS{T}(br::AbstractAlgebra.Ring, a::RingElement) where T <: RingElement
        return new(br, br(a))
    end

    function ConstPS{T}(a::T) where T <: RingElement
        return new(parent(a), a)
    end
end

function get_truncation(ps::ConstPS{T}, R::Generic.SeriesRing{T}) where T <: RingElement
    return R(ps.c)
end

#------------------------------------------------------------------------------

mutable struct RationalPS{T} <: ComputablePS{T}
    base_ring::AbstractAlgebra.Ring
    numer::PolyElem{T}
    denom::PolyElem{T}
    cached::Array{T, 1}

    function RationalPS{T}(numer::PolyElem{T}, denom::PolyElem{T}) where T <: RingElement
        return new(base_ring(parent(numer)), numer, denom, Array{T, 1}())
    end
end

function get_truncation_raw(ps::RationalPS, R::Generic.SeriesRing)
    t = gen(R)
    return evaluate(ps.numer, t) * inv(evaluate(ps.denom, t))
end

#------------------------------------------------------------------------------

struct IdentityPS{T} <: ComputablePS{T}
    base_ring::AbstractAlgebra.Ring
    cached::Array{<: T, 1}

    function IdentityPS{T}(br::AbstractAlgebra.Ring) where T <: RingElement
        return new(br, Array{T, 1}())
    end
end

function get_truncation(ps::IdentityPS, R::Generic.SeriesRing)
    return gen(R)
end

#------------------------------------------------------------------------------

# log(1 + t)
mutable struct LogPS{T} <: ComputablePS{T}
    base_ring::AbstractAlgebra.Ring
    cached::Array{<: T, 1}

    function LogPS{T}(br::AbstractAlgebra.Ring) where T <: RingElement
        return new(br, Array{T, 1}())
    end
end

function get_truncation_raw(ps::LogPS, R::Generic.SeriesRing)
    t = gen(R)
    return sum([(-1)^(i - 1) // i * t^i for i in 1:max_precision(R)])
end

#------------------------------------------------------------------------------

mutable struct ExpPS{T} <: ComputablePS{T}
    base_ring::AbstractAlgebra.Ring
    cached::Array{<: T, 1}

    function ExpPS{T}(br::AbstractAlgebra.Ring) where T <: RingElement
        return new(br, Array{T, 1}())
    end
end

function get_truncation_raw(ps::ExpPS, R::Generic.SeriesRing)
    t = gen(R)
    return sum([1 // factorial(big(i)) * t^i for i in 0:max_precision(R)])
end

#------------------------------------------------------------------------------
# Arithmetic of computable power series
#------------------------------------------------------------------------------

mutable struct CompositePS{T} <: ComputablePS{T}
    base_ring::AbstractAlgebra.Ring
    func # function to be applied
    left::ComputablePS{T}
    right::ComputablePS{T}
    cached::Array{<: T, 1}

    function CompositePS{T}(func, left, right) where T <: RingElement
        return new(base_ring(left), func, left, right, Array{T, 1}())
    end
end

function get_truncation_raw(ps::CompositePS, R::Generic.SeriesRing)
    lps = get_truncation(ps.left, R)
    rps = get_truncation(ps.right, R)
    return ps.func(lps, rps)
end

#------------------------------------------------------------------------------

mutable struct CompositeUnaryPS{T} <: ComputablePS{T}
    base_ring::AbstractAlgebra.Ring
    func
    arg::ComputablePS{T}
    cached::Array{<: T, 1}

    function CompositeUnaryPS{T}(func, arg::ComputablePS{T}) where T <: RingElement
        return new(base_ring(arg), func, arg, Array{T, 1}())
    end
end

function get_truncation_raw(ps::CompositeUnaryPS, R::Generic.SeriesRing)
    inner_ps = get_truncation(ps.arg, R)
    return ps.func(inner_ps)
end

#------------------------------------------------------------------------------

function Base.:-(f::ComputablePS{T}) where T <: RingElement
    return CompositeUnaryPS{T}(a -> -a, f)
end

function Base.:^(f::ComputablePS{T}, n::Integer) where T <: RingElement
    return CompositeUnaryPS{T}(a -> a^n, f)
end

#------------------------------------------------------------------------------

function (f::ComputablePS{T})(g::ComputablePS{T}) where T <: RingElement
    return CompositePS{T}((a, b) -> ps_compose(a, b), f, g)
end

function Base.:+(f::ComputablePS{T}, g::ComputablePS{T}) where T <: RingElement
    return CompositePS{T}((a, b) -> a + b, f, g)
end

function Base.:*(f::ComputablePS{T}, g::ComputablePS{T}) where T <: RingElement
    return CompositePS{T}((a, b) -> a * b, f, g)
end

function Base.:-(f::ComputablePS{T}, g::ComputablePS{T}) where T <: RingElement
    return CompositePS{T}((a, b) -> a - b, f, g)
end

function Base.:/(f::ComputablePS{T}, g::ComputablePS{T}) where T <: RingElement
    return CompositePS{T}((a, b) -> a * inv(b), f, g)
end

function AbstractAlgebra.parent(ps::ComputablePS{T}) where T <: RingElement
    return ConstPS{T}
end

AbstractAlgebra.promote_rule(::Type{ComputablePS{T}}, ::Type{T}) where T <: RingElement = ComputablePS{T}
AbstractAlgebra.promote_rule(::Type{T}, ::Type{ComputablePS{T}}) where T <: RingElement = ComputablePS{T}


#------------------------------------------------------------------------------
# Defining shifts and other properties for nondifference rings and constants
#------------------------------------------------------------------------------

function sigma(a::RingElement)
    return a
end

function sigma(a::RingElement, n::Integer)
    if n == 0
        return a
    end
    return sigma(sigma(a, n - 1))
end

function get_truncation(a::RingElement, R::Generic.SeriesRing)
    return R(a)
end

#------------------------------------------------------------------------------
# Difference polynomials and operations on them
#------------------------------------------------------------------------------
# All polynomials here are assumed to be polynomials in x0, ..., xMAX_ORD
# over a SigmaPSDomain or a ring of constants

function sigma(poly::MPolyElem)
    builder = MPolyBuildCtx(parent(poly))
    for term in zip(exponent_vectors(poly), coeffs(poly))
        exp, coef = term
        if exp[end] != 0
            throw(Base.ArgumentError("You are exceeding the maximal order $MAX_ORD, increase it and try again"))
        end
        pushfirst!(exp, 0)
        pop!(exp)
        push_term!(builder, sigma(coef), exp)
    end
    res = finish(builder)
    return res
end

function sigma(poly::MPolyElem, i::Integer)
    if i == 0
        return poly
    end
    return sigma(sigma(poly, i - 1))
end

# used for reduce! of SigmaPSDomainElem
function simplify(poly::MPolyElem)
    builder = MPolyBuildCtx(parent(poly))
    for term in zip(exponent_vectors(poly), coeffs(poly))
        exp, coef = term
        if !iszero(coef)
            push_term!(builder, coef, exp)
        end
    end
    return finish(builder)
end

#------------------------------------------------------------------------------

function order(poly::MPolyElem)
    result = 0
    for exp in exponent_vectors(poly)
        if sum(exp) != 0
            result = max(result, findlast(map(x -> x > 0, exp)))
        end
    end
    return result - 1
end

function rittrank(poly::MPolyElem)
    ord = order(poly)
    if ord == -1
        return (-1, -1)
    end
    return (ord, degree(poly, ord + 1))
end

function leader(poly::MPolyElem)
    ord = order(poly)
    if ord == -1
        throw(ArgumentError("Cannot compute the leader of a constant $poly"))
    end
    return gens(parent(poly))[ord + 1]
end

#------------------------------------------------------------------------------

function initial(poly::MPolyElem)
    ord = order(poly)
    if ord == -1
        throw(ArgumentError("Cannot compute the initial of a constant"))
    end
    d = degree(poly, gens(parent(poly))[ord + 1])
    
    builder = MPolyBuildCtx(parent(poly))
    for term in zip(exponent_vectors(poly), coeffs(poly))
        exp, coef = term
        if exp[ord + 1] == d
            new_exp = [(i == ord + 1 ? 0 : exp[i]) for i in 1:length(exp)]
            push_term!(builder, coef, new_exp)
        end
    end
    res = finish(builder)
    return res
end

function get_truncation(a::fmpq, n::Integer)
    return a
end

function separant(poly::MPolyElem)
    # TODO: use derivative function
    ind = order(poly) + 1
    builder = MPolyBuildCtx(parent(poly))
    for term in zip(exponent_vectors(poly), coeffs(poly))
        exp, coef = term
        if exp[ind] > 0
            new_coef = coef * exp[ind]
            new_exp = [(i == ind ? exp[i] - 1 : exp[i]) for i in 1:length(exp)]
            push_term!(builder, new_coef, new_exp)
        end
    end
    res = finish(builder)
    return res
end

#------------------------------------------------------------------------------

function isreducible(p::MPolyElem, q::MPolyElem)
    prank = rittrank(p)
    qrank = rittrank(q)
    return (prank[1] >= qrank[1]) && (prank[2] >= qrank[2])
end

function isreducible(p::MPolyElem, Q::Array{<: MPolyElem, 1})
    return any([isreducible(p, q) for q in Q])
end

function reduce(p::MPolyElem, q::MPolyElem)
    if !isreducible(p, q)
        return p
    end
    shift_count = order(p) - order(q)
    qs = sigma(q, shift_count)

    l = leader(p)
    Iqs = initial(qs)
    while degree(p, l) >= degree(qs, l)
        p = Iqs * p - initial(p) * l^(degree(p, l) - degree(qs, l)) * qs
    end
    return reduce(p, q)
end

function reduce(p::MPolyElem, chain::Array{<: MPolyElem, 1})
    if length(chain) == 0
        return p
    end
    return reduce(reduce(p, chain[end]), chain[1:end - 1])
end

#------------------------------------------------------------------------------

function delta_polynomial(p::MPolyElem, q::MPolyElem)
    if order(p) > order(q)
        return delta_polynomial(q, p)
    end
    diff_ord = order(q) - order(p)
    diff_deg = degree(p, leader(p)) - degree(q, leader(q))
    return initial(q) * sigma(p, diff_ord) - sigma(initial(p), diff_ord) * leader(q)^diff_deg * q
end

function select_charset(polys::Array{<: MPolyElem, 1})
    withranks = map(p -> (rittrank(p), p), polys)
    sort!(withranks, by = x -> x[1])
    result = Any[withranks[1]]
    for x in withranks[2:end]
        rank = x[1]
        prev = result[end][1]
        if !((rank[1] >= prev[1]) && (rank[2] >= prev[2]))
            push!(result, x)
        end
    end
    return map(x -> x[2], result)
end

#------------------------------------------------------------------------------
# Difference polynomials: evaluations, valuations, and indicial polynomials
#------------------------------------------------------------------------------

SIGMA_TO_DELTA = Dict()

function sigma_to_delta(shift::ComputablePS{T}) where T <: RingElement
    if length(SIGMA_TO_DELTA) != 0
        return SIGMA_TO_DELTA[1]
    end
    k =  valuation(shift - IdentityPS{T}(base_ring(shift)))
    denom = IdentityPS{T}(base_ring(shift))^(k - 1)
    result = []
    push!(result, ComputablePS[(i == 0 ? ConstPS{T}(base_ring(shift), 1) : ConstPS{T}(base_ring(shift), 0)) for i in 0:MAX_ORD])
    for i in 1:MAX_ORD
        prev = result[end]
        next = copy(prev)
        for j in i:-1:1
            next[j + 1] = next[j + 1] + denom * next[j](shift)
            next[j] = next[j](shift)
        end
        push!(result, next)
    end
    SIGMA_TO_DELTA[1] = result
    return result
end

function diff_evaluate(poly::MPolyElem, f::ComputablePS{T}, shift::ComputablePS{T}) where T <: RingElement
    max_ord = length(gens(parent(poly))) - 1
    result = ConstPS{T}(base_ring(f), 0)
    gen_shifts = ComputablePS{T}[f]
    for i in 1:max_ord
        push!(gen_shifts, gen_shifts[end](shift))
    end

    for term in zip(exponent_vectors(poly), coeffs(poly))
        exp, coef = term
        result = result + coef * prod([gen_shifts[i]^exp[i] for i in 1:(max_ord + 1)])
    end
    return result
end

function fvaluation(poly::MPolyElem, f::ComputablePS, shift::ComputablePS)
    eval = diff_evaluate(poly, f, shift)
    return valuation(eval)
end

function sigma_linear_part(poly::MPolyElem, f::ComputablePS, shift::ComputablePS)
    ord = order(poly)
    result = []
    for i in 0:ord
        push!(result, diff_evaluate(derivative(poly, gens(parent(poly))[i + 1]), f, shift))
    end
    return result
end

function delta_linear_part(poly::MPolyElem, f::ComputablePS{T}, shift::ComputablePS{T}) where T <: RingElement
    std = sigma_to_delta(shift)
    slp = sigma_linear_part(poly, f, shift)

    l = length(slp)
    result = ComputablePS{T}[ConstPS{T}(base_ring(f), 0) for i in 1:l]
    for i in 1:l
        for j in 1:l
            result[j] = result[j] + slp[i] * std[i][j]
        end
    end
    return result
end

function valuation_linear_part(poly::MPolyElem, f::ComputablePS, shift::ComputablePS)
    dlp = delta_linear_part(poly, f, shift)
    return valuation(dlp)[1]
end

function indpoly_linear_part(poly::MPolyElem, f::ComputablePS{T}, shift::ComputablePS{T}) where T <: RingElement
    dlp = delta_linear_part(poly, f, shift)
    val = valuation_linear_part(poly, f, shift)
    prec = val + 1
    R, t = PowerSeriesRing(base_ring(f), prec, "t"; model=:capped_absolute)
    eval_dlp = [get_truncation(ps, R) for ps in dlp]
    coeffs = [coeff(ps, val) for ps in eval_dlp]

    k =  valuation(shift - IdentityPS{T}(base_ring(f)))
    R, t = PowerSeriesRing(base_ring(f), k + 1, "t"; model=:capped_absolute)
    g = coeff(get_truncation(shift, R), k)
    
    R, x = PolynomialRing(base_ring(f), "x")
    return sum([coeffs[i] * (g * x)^(i - 1) for i in 1:length(coeffs)])
end

function Z_linear_part(poly::MPolyElem, f::ComputablePS, shift::ComputablePS)
    ind_poly = indpoly_linear_part(poly, f, shift)
    return max(-1, integer_roots(ind_poly)...)
end

#------------------------------------------------------------------------------
# Computable power series defined by an annihilator and IC
#------------------------------------------------------------------------------

mutable struct AnnihilatorBasedPS{T} <: ComputablePS{T}
    base_ring::AbstractAlgebra.Ring
    annihilator::MPolyElem
    shift::ComputablePS{T}
    known_terms::Array{Any, 1}

    function AnnihilatorBasedPS{T}(annihilator, shift, known_terms) where T <: RingElement
        return new(base_ring(shift), annihilator, shift, known_terms)
    end
end

function get_truncation(ps::AnnihilatorBasedPS{T}, R::Generic.SeriesRing{T}) where T <: RingElement
    S, x = PolynomialRing(base_ring(ps), "x")
    trunc = sum([ps.known_terms[i] * x^(i - 1) for i in 1:length(ps.known_terms)])
    prec = max_precision(R)
    while length(ps.known_terms) < prec
        f = RationalPS{T}(trunc, S(1))
        mu = valuation_linear_part(ps.annihilator, f, ps.shift)
        ind_poly = indpoly_linear_part(ps.annihilator, f, ps.shift)
        n = length(ps.known_terms)
        Pf = get_truncation(diff_evaluate(ps.annihilator, f, ps.shift), n + mu + 1)
        push!(ps.known_terms, -coeff(Pf, n + mu) // evaluate(ind_poly, n))
        trunc += ps.known_terms[end] * x^n
    end
    t = gen(R)
    return sum([ps.known_terms[i + 1] * t^i for i in 0:(prec - 1)])
end

#------------------------------------------------------------------------------
# Sigma domains
#------------------------------------------------------------------------------

struct SigmaPSDomain{T} <: AbstractSigmaPSDomain
    coef_ring::AbstractAlgebra.Ring
    shift::ComputablePS{T}
    base_ring::AbstractAlgebra.Ring
    generator::ComputablePS{T}
    diff_poly_ring::AbstractAlgebra.MPolyRing
    annihilator::AbstractAlgebra.MPolyElem
    sigma_depth::Integer
end

function AbstractAlgebra.base_ring(R::SigmaPSDomain)
    return R.base_ring
end

function sigma_depth(R::SigmaPSDomain)
    return R.sigma_depth
end

#------------------------------------------------------------------------------

function construct_domain(shift::ComputablePS{T}, ps::ComputablePS{T}) where T <: RingElement
    cring = base_ring(ps)
    diff_poly_ring, v = PolynomialRing(cring, ["x$i" for i in 0:MAX_ORD])
    return SigmaPSDomain{T}(cring, shift, cring, ps, diff_poly_ring, zero(diff_poly_ring), 0)
end

function adjoin_transcendental(base::SigmaPSDomain{T}, ps::ComputablePS{<: T}) where T <: RingElement
    diff_poly_ring, v = PolynomialRing(base, ["x$i" for i in 0:MAX_ORD])
    return SigmaPSDomain{T}(base.coef_ring, base.shift, base, ps, diff_poly_ring, zero(diff_poly_ring), sigma_depth(base) + 1)
end

#------------------------------------------------------------------------------

function construct_domain(ps::AnnihilatorBasedPS{T}) where T <: RingElement
    cring = base_ring(ps)
    diff_poly_ring, v = PolynomialRing(cring, ["x$i" for i in 0:MAX_ORD])
    annihilator = parent_ring_change(ps.annihilator, diff_poly_ring)
    return SigmaPSDomain{T}(cring, ps.shift, cring, ps, diff_poly_ring, annihilator, 0)
end

function adjoin_sigma_algebraic(base::SigmaPSDomain{T}, ps::AnnihilatorBasedPS{<: T}) where T <: RingElement
    diff_poly_ring, v = PolynomialRing(base, ["x$i" for i in 0:MAX_ORD])
    annihilator = parent_ring_change(ps.annihilator, diff_poly_ring)
    return SigmaPSDomain{T}(base.coef_ring, base.shift, base, ps, diff_poly_ring, annihilator, sigma_depth(base) + 1)
end

#------------------------------------------------------------------------------
# Sigma domain elements
#------------------------------------------------------------------------------

mutable struct SigmaPSDomainElem{T} <: ComputablePS{T}
    base_ring::AbstractAlgebra.Ring
    parent::SigmaPSDomain{T}
    poly::MPolyElem

    function SigmaPSDomainElem{T}(parent::SigmaPSDomain{T}, poly::MPolyElem) where T <: RingElem
        return new(parent.coef_ring, parent, poly)
    end
end

AbstractAlgebra.parent_type(::Type{SigmaPSDomainElem{T}}) where T <: RingElem = SigmaPSDomain{T}
AbstractAlgebra.elem_type(::Type{SigmaPSDomain{T}}) where T <: RingElem = SigmaPSDomainElem{T}

function AbstractAlgebra.parent(a::SigmaPSDomainElem)
    return a.parent
end

function Base.show(ps::SigmaPSDomainElem)
    show(poly(ps))
end

function AbstractAlgebra.needs_parentheses(ps::SigmaPSDomainElem)
    return true
end

function poly(a::SigmaPSDomainElem)
    return a.poly
end

AbstractAlgebra.promote_rule(::Type{T}, ::Type{SigmaPSDomainElem{T}}) where T <: RingElement = SigmaPSDomainElem{T}
AbstractAlgebra.promote_rule(::Type{SigmaPSDomainElem{T}}, ::Type{T}) where T <: RingElement = SigmaPSDomainElem{T}
AbstractAlgebra.promote_rule(::Type{Generic.MPolyElem}, ::Type{SigmaPSDomainElem}) = Generic.MPolyElem
AbstractAlgebra.promote_rule(::Type{SigmaPSDomainElem}, ::Type{Generic.MPolyElem}) = Generic.MPolyElem

function (a::AbstractAlgebra.Generic.MPolyRing{SigmaPSDomainElem{T}})(b::SigmaPSDomainElem{T}) where {T <: RingElement}
    b = base_ring(a)(b)
    z = AbstractAlgebra.Generic.MPoly{SigmaPSDomainElem{T}}(a, b)
    return z
end

#------------------------------------------------------------------------------

function Base.zero(R::SigmaPSDomain{T}) where T <: RingElem
    return SigmaPSDomainElem{T}(R, zero(R.diff_poly_ring))
end

function Base.one(R::SigmaPSDomain{T}) where T <: RingElem
    return SigmaPSDomainElem{T}(R, one(R.diff_poly_ring))
end

function (R::SigmaPSDomain)()
    return zero(R)
end

#------------------------------------------------------------------------------

function AbstractAlgebra.gen(R::SigmaPSDomain{T}) where T <: RingElem
    return SigmaPSDomainElem{T}(R, gens(R.diff_poly_ring)[1])
end

function (R::SigmaPSDomain{T})(a::SigmaPSDomainElem) where T <: RingElement
    if parent(a) == R
        return a
    end
    if parent(a) == base_ring(R)
        return SigmaPSDomainElem{T}(R, R.diff_poly_ring(a))
    end
    return R(base_ring(R)(a))
end

function (R::SigmaPSDomain{T})(a::RingElement) where T <: RingElement
    if base_ring(R) == R.coef_ring
        return SigmaPSDomainElem{T}(R, R.diff_poly_ring(a))
    end
    return R(base_ring(R)(a))
end

#------------------------------------------------------------------------------

function common_parent(a::SigmaPSDomainElem, b::SigmaPSDomainElem)
    if sigma_depth(parent(a)) > sigma_depth(parent(b))
        return parent(a)
    end
    return parent(b)
end

function Base.:+(a::SigmaPSDomainElem{T}, b::SigmaPSDomainElem{T}) where T <: RingElement
    p = common_parent(a, b)
    return SigmaPSDomainElem{T}(p, poly(p(a)) + poly(p(b)))
end

function Base.:*(a::SigmaPSDomainElem{T}, b::SigmaPSDomainElem{T}) where T <: RingElement
    p = common_parent(a, b)
    return SigmaPSDomainElem{T}(p, poly(p(a)) * poly(p(b)))
end

function Base.:-(a::SigmaPSDomainElem{T}, b::SigmaPSDomainElem{T}) where T <: RingElement
    p = common_parent(a, b)
    return SigmaPSDomainElem{T}(p, poly(p(a)) - poly(p(b)))
end

function Base.:^(a::SigmaPSDomainElem{T}, n::Integer) where T <: RingElement
    return SigmaPSDomainElem{T}(parent(a), poly(a)^n)
end

function Base.:-(a::SigmaPSDomainElem{T}) where T <: RingElement
    return SigmaPSDomainElem{T}(parent(a), -poly(a))
end

function Base.:(==)(a::SigmaPSDomainElem{T}, b::SigmaPSDomainElem{T}) where T <: RingElement
    p = common_parent(a, b)
    return poly(p(a)) == poly(p(b))
end

function Base.:(!=)(f::SigmaPSDomainElem, g::SigmaPSDomainElem)
    return !iszero(f - g)
end

function AbstractAlgebra.iszero(p::MPolyElem{SigmaPSDomainElem{T}}) where T <: RingElem
    for c in coeffs(p)
        if !iszero(c)
            return false
        end
    end
    return true
end

function AbstractAlgebra.reduce!(p::SigmaPSDomainElem{T}) where T <: RingElem
    p.poly = simplify(p.poly)
    return p
end
#------------------------------------------------------------------------------

function get_truncation(ps::SigmaPSDomainElem, R::Generic.SeriesRing)
    return get_truncation(diff_evaluate(poly(ps), parent(ps).generator, parent(ps).shift), R)
end

function sigma(ps::SigmaPSDomainElem{T}) where T <: RingElem
    return SigmaPSDomainElem{T}(parent(ps), sigma(poly(ps)))
end

#------------------------------------------------------------------------------
# Zero testing
# -----------------------------------------------------------------------------

FAST_CHECK_ORDER = 5

function zerotest_inner(f::SigmaPSDomainElem{T}, Q::Array{<: MPolyElem, 1}) where T <: RingElement
    Q = map(simplify, Q)
    if any(map(p -> total_degree(p) == 0 && !iszero(first(coeffs(p))), Q))
        return false
    end
    Q = filter(p -> total_degree(p) > 0, Q)

    if length(Q) == 0
        return true
    end

    # fast initial check
    if FAST_CHECK_ORDER > 0
        check_ord = (sigma_depth(parent(f)) == 0 ? FAST_CHECK_ORDER * 2 : FAST_CHECK_ORDER)
        val = map(p -> AbstractAlgebra.valuation(get_truncation(diff_evaluate(p, f, parent(f).shift), check_ord)), Q)
        if min(val...) < check_ord
            return false
        end
    end

    charset = select_charset(Q)

    for r in charset
        for si in [initial(r), separant(r)]
            if total_degree(si) == 0
                continue
            end
            rem = reduce(si, charset)
            if !iszero(rem)
                res = zerotest_inner(f, [si, Q...])
                if res
                    return true
                else
                    val = valuation(map(p -> diff_evaluate(p, f, parent(f).shift), [si, Q...]))
                    if val[2] != 1
                        return false
                    end
                end
            end
        end
    end

    for q in [parent(f).annihilator, Q...]
        rem = reduce(q, charset)
        if !iszero(rem)
            return zerotest_inner(f, [rem, Q...])
        end
    end

    for i in 2:length(charset)
        dpoly = delta_polynomial(charset[i - 1], charset[i])
        rem = reduce(dpoly, charset)
        if !iszero(rem)
            return zerotest_inner(f, [rem, Q...])
        end
    end

    s = 0
    for i in 2:length(charset)
        cur = charset[i]
        degree_diff = degree(charset[i - 1], leader(charset[i - 1])) - degree(cur, leader(cur)) + 1
        val_init = fvaluation(initial(cur), f, parent(f).shift)
        val_sep = fvaluation(separant(cur), f, parent(f).shift)
        s = max(
            s,
            degree_diff * val_init + val_sep
        )
    end

    s = max(
        s,
        valuation_linear_part(parent(f).annihilator, f, parent(f).shift),
        Z_linear_part(parent(f).annihilator, f, parent(f).shift),
        valuation_linear_part(charset[1], f, parent(f).shift)
    )

    prec = Int(numerator(2 * s + 1))
    evals = [get_truncation(diff_evaluate(r, f, parent(f).shift), prec) for r in charset]
    return !any(map(ps -> AbstractAlgebra.valuation(ps) < prec, evals))
end

function AbstractAlgebra.iszero(f::SigmaPSDomainElem)
    if iszero(poly(f))
        return true
    end
    if iszero(parent(f).annihilator)
        return false
    end
    if (total_degree(poly(f)) == 0 && iszero(first(coeffs(poly(f))))) || zerotest_inner(gen(parent(f)), [poly(f)])
        return true
    end
    return false
end
