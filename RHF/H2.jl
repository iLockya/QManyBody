include("gaussian.jl")
include("RHFSCF.jl")
using LinearAlgebra
using TensorOperations
using ArgParse

s = ArgParseSettings()
@add_arg_table! s begin
    "-d"
        help = "Two nuclear distance."
        arg_type = Float64
        default = 1.4
    "--iter", "-i"
        help = "max iteration."
        arg_type = Int64
        default = 100
    "--eps", "-e"
        help = "convergence criterial."
        arg_type = Float64
        default = 1e-9
    "--criteria", "-c"
        help = "convergence criterial.'energy' or 'density'. "
        arg_type = String
        default = "energy"
    "--prt", "-p"
        help = "print the output of SCF procedure."
        action = :store_true
end
args = parse_args(s)



### Run RHF SCF.
system = init_H2(args["d"])
RHFSCF!(system; max_iter=args["iter"],eps=args["eps"],criterion=args["criteria"],prt=args["prt"])


### Calculate H energy by STO-3G basis, -0.5 for theoretical value.
E₀ = ∫(STO3G(),STO3G();type = "kinetic") + ∫(STO3G(),STO3G();type = "nucleus") # -0.46658185041114275
Bond_Energy = system.E - 2E₀

### print information.
println("Orbital 1: |ψ₁⟩ = $(system.Orbitals[1,1])|φ₁⟩ + $(system.Orbitals[2,1])|φ₂⟩, energy = $(system.ε[1]).")
println("Orbital 2: |ψ₂⟩ = $(system.Orbitals[1,2])|φ₁⟩ + $(system.Orbitals[2,2])|φ₂⟩, energy = $(system.ε[2]).")
println("Total Energy at $(args["d"]) a.u. is $(system.E).")
println("Bond Energy at $(args["d"]) a.u. is $Bond_Energy.")