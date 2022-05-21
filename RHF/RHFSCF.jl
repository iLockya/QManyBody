using LinearAlgebra
using TensorOperations

include("gaussian.jl")

mutable struct Molecule{T,S}
    N::T # Number of electrons. N is even in Closed-Shell.
    Basis::Vector{<:STO} # Basis set.
    Nus::Vector{Nucleus{T,S}} # Nucleus imformation.
    Orbitals::Matrix{S} # each colum is an orbital coefficients.
    ε::Vector{S} # Energy level of each orbital.
    E::S # Total energy, include electrons energy and nuclear energy.
end

function init_H2(d::Float64; Orbitals = zeros(2,2))
    """
    Construct H₂.
    """
    Basis = [STO3G(R=[0,0.0,0.0]),STO3G(R=[d,0.0,0.0])]
    Nus = [Nucleus(1,[0.0,0.0,0.0]),Nucleus(1,[d,0.0,0.0])]
    ε = zeros(2)
    system = Molecule(2,Basis,Nus,Orbitals,ε,Inf)
    return system
end

function Basis_Overlap_Matrix(Basis;type="overlap",N::Nucleus=Nucleus())
    n = length(Basis)
    S = zeros(n,n)
    for i=1:n
        for j =i:n
            if i==j
                S[i,i] = ∫(Basis[i],Basis[i];type,N)
            else
                S[i,j] = ∫(Basis[i],Basis[j];type,N)
                S[j,i] = S[i,j]
            end
        end
    end
    return S
end

function H_core(Basis,Nus)
    """
    Single particle Hamiltonian.
    """
    T = Basis_Overlap_Matrix(Basis;type="kinetic")
    V = [Basis_Overlap_Matrix(Basis;type="nucleus",N=N) for N ∈ Nus] |> sum
    return T + V
end

function e_eInteration_Matrix(Basis)
    n = length(Basis)
    g = zeros(n,n,n,n)
    for μ=1:n,ν=1:n,σ=1:n,λ=1:n
        g[μ,ν,σ,λ] = ∫(Basis[μ],Basis[ν],Basis[σ],Basis[λ]) - 
                     0.5* ∫(Basis[μ],Basis[λ],Basis[σ],Basis[ν])
    end
    return g
end

function n_nEnergy(Nus)
    """
    Nuclear repulsion energy.
    """
    n = length(Nus)
    E = 0
    for i=1:n,j=i+1:n
        E += Nus[i].Z*Nus[j].Z / norm(Nus[i].R - Nus[j].R)
    end
    return E
end

# function S_normalize!(C,S)
#     for i=1:size(C,2)
#         N = C[:,i]'*S*C[:,i]
#         C[:,i] = C[:,i] / sqrt(N)
#     end
# end

function RHFSCF!(system::Molecule; max_iter::Int64 = 100,eps::Float64 = 1e-9, 
                              criterion::String = "energy",prt::Bool = false)
    """
    Run RHF SCF procedure.
    system: System object.
    mat_iter: maximum number of iterations, default is 100.
    eps: convergence criteria, default is 1e-6.
    criterion: "energy" for the lowest orbital energy or "density" for the 
               density matrix, default is "energy".
    prt: print the output of SCF procedure, default is false.
    """

    ### Step 1. Calculate Overlap Matrix S, single e⁻ Hamiltonian h and e-e interaction tensor g.
    ###         These values only need to calulate once.

    # Overlap Matrix.
    S = Basis_Overlap_Matrix(system.Basis) 
    s,U = eigen(S)
    X = U / Diagonal(sqrt.(s))

    # Single e⁻ Hamiltonian
    h = H_core(system.Basis,system.Nus)

    # e-e interaction tensor
    g = e_eInteration_Matrix(system.Basis) 

    # Nuclear repulsion energy.
    En_n = n_nEnergy(system.Nus)

    for i=1:max_iter

        ### Step 2. Calculate density matrix D and e-e interaction matrix.
        D = [2*system.Orbitals[:,j]*system.Orbitals[:,j]' for j in 1:system.N÷2] |> sum

        ### Step 3. Calculate Fock matrix.
        @tensor F[μ,ν] := h[μ,ν] + g[μ,ν,σ,λ]*D[σ,λ]

        ### Step 4. Diagnolize F under the overlap matrix S.
        # system.ε, system.Orbitals = eigen(F,S)
        # system.ε, system.Orbitals = eigen(inv(S)*F)
        # S_normalize!(system.Orbitals,S)
        system.ε, system.Orbitals = eigen(X'*F*X)
        system.Orbitals = X * system.Orbitals

        ### Step 5. Check Convergence criterion.
        if criterion == "energy"
            E_ = sum( D.*(h+F)) / 2 + En_n
            δ = abs(system.E-E_)
        else
            D_ = [2*system.Orbitals[:,j]*system.Orbitals[:,j]' for j in 1:system.N÷2] |> sum
            δ = norm(D - D_) / length(system.Basis)
        end

        ## Update total electrons energy. 
        system.E = sum(D.*(h+F))/2 + En_n
        if prt
            println("#Step $i:")
            println("   Orbital 1: |ψ₁⟩ = $(system.Orbitals[1,1])|φ₁⟩ + $(system.Orbitals[2,1])|φ₂⟩, energy = $(system.ε[1]).")
            println("   Orbital 2: |ψ₂⟩ = $(system.Orbitals[1,2])|φ₁⟩ + $(system.Orbitals[2,2])|φ₂⟩, energy = $(system.ε[2]).")
            println("   Total energy is $(system.E).")
        end

        if δ < eps
            if prt
                println("Converged at step $(i)!")
            end
            break
        end
    end
end


