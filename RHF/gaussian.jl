using SpecialFunctions
import Base.*


struct Gaussian{T}
    """
    Gaussian type function.
    χ(r) = (2α/π)^(3/4) * e^(-α|r-R|^2)
    The Gaussian type function is normalized in this code.
    """
    α::T
    R::Vector{T}
end

@Base.kwdef struct Nucleus{T,S}
    Z::T = 1
    R::Vector{S} = [0.0,0.0,0.0]
end

abstract type STO end
@Base.kwdef struct STO3G{T} <: STO
    """
    STO-3G basis set for H 1s by default.
    """
    n=3
    d::Vector{T} = [0.4446345422,0.5353281423,0.1543289673]
    α::Vector{T} = [0.1688554040,0.6239137298,3.425250914]
    R::Vector{T} = [0.0,0.0,0.0]
    #φ::Vector{Gaussian{T}} = [Gaussian(0.168856,[0.0,0.0,0.0]),
    #                          Gaussian(0.6213913,[0.0,0.0,0.0]),
    #                          Gaussian(3.42525,[0.0,0.0,0.0])]
end

@Base.kwdef struct STO6G{T} <: STO
    """
    STO-6G basis set for H 1s by default.
    """
    n=6
    d::Vector{T} = [0.009163596281,0.04936149294,0.1685383049,0.3705627997,0.4164915298,0.1303340841]
    α::Vector{T} = [35.52322122,6.513143725,1.822142904,0.6259552659,0.2430767471,0.1001124280]
    R::Vector{T} = [0.0,0.0,0.0]
    #φ::Vector{Gaussian{T}} = [Gaussian(0.168856,[0.0,0.0,0.0]),
    #                          Gaussian(0.6213913,[0.0,0.0,0.0]),
    #                          Gaussian(3.42525,[0.0,0.0,0.0])]
end


function *(A::Gaussian,B::Gaussian)
    """
    product of two Gaussians.
    return a Gaussian.
    """
    p = A.α + B.α
    R = (A.α.*A.R + B.α.*B.R)/p
    K = ( 2*A.α*B.α/(p*π) )^(3/4) * exp(-A.α*B.α*sum( (A.R-B.R).^2 )/p)
    return K,Gaussian(p,R)
end


function ∫(A::Gaussian,B::Gaussian; type::String = "overlap",N::Nucleus=Nucleus())  
    """
    two center integral.
    overlap:  ⟨A|B⟩
    kinetic:  ⟨A|-1/2∇²|B⟩.
    potential: ⟨A|1/r|B⟩
    """
    @assert type ∈ ["overlap","kinetic","potential","nucleus"]
    p = A.α + B.α
    q = A.α*B.α*sum( (A.R-B.R).^2 )/p
    if type == "overlap"
        return ( 4*A.α*B.α/(p*p) )^(3/4) * exp(-q) 
    elseif type == "kinetic"
        return 2^(3/2)*(A.α*B.α)^(7/4) / p^(5/2) * (3-2q) * exp(-q)
    elseif type == "potential"
        return 2^(5/2)*π/ ( (A.α*B.α)^(1/4) * p^(1/2)) * F₀(q)
    else
        Rp = (A.α.*A.R + B.α.*B.R)/p
        t = p*sum( (Rp-N.R).^2 )
        return -(A.α*B.α)^(3/4)*2^(5/2)*N.Z/(p*π^(1/2)) * exp(-q) * F₀(t)
    end
end

function ∫(A::STO,B::STO;type::String = "overlap",N::Nucleus=Nucleus())
    s = 0
    for i = 1:A.n
        for j = i:B.n
            if j == i
                s += A.d[i]*B.d[i] * ∫(Gaussian(A.α[i],A.R),Gaussian(B.α[i],B.R);type,N)
            else
                s += 2*A.d[i]*B.d[j] * ∫(Gaussian(A.α[i],A.R),Gaussian(B.α[j],B.R);type,N)
            end
        end
    end
    return s
end


function ∫(A::Gaussian,B::Gaussian,C::Gaussian,D::Gaussian)
    """
    e-e interaction.
    """
    ## First, combine A&B, C&D.
    K₁, AB = A*B
    K₂, CD = C*D
    ## Second, conver to two Gaussian integral.
    return K₁*K₂*∫(AB,CD;type="potential")
end

function ∫(A::STO,B::STO,C::STO,D::STO)
    """
    e-e interaction.
    """
    s=0
    for i=1:A.n,j=1:B.n,k=1:C.n,l=1:D.n
        s += A.d[i]*B.d[j]*C.d[k]*D.d[l] * ∫(Gaussian(A.α[i],A.R),Gaussian(B.α[j],B.R),Gaussian(C.α[k],C.R),Gaussian(D.α[l],D.R))
    end
    return s
end


function F₀(t::Float64)
    if t < 1e-9
        return 1
    else
        return 1/2*(π/t)^(1/2)*erf(t^(1/2))
    end
end


