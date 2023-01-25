include("functions.jl")
using BoundaryValueDiffEq
using DifferentialEquations
using Printf
using Plots

function modelo_marc!(du, u, p, t)
	θ = u[1]
	dθ = u[2]
	@printf "Guess: %.4f °C em %.4f mm\n" θ t
	du[1] = dθ
	du[2] = 1/𝑘(θ) * (-w * ρ(θ) * cv(θ) * dθ + 𝒉(t) * (θ - θ∞) - j^2 * r)
end
function bc!(residual, u, p, t)
	residual[1] = u[1][1] - 1515.
	residual[2] = u[end][1] - 158.
end

const w = 0.075; const j = 99.47; const ΔH = 120.; const r = 5.62e-6 
const Va = 20; const θ∞ = 25
contorno = ((0.,1515.),(15.,158.))
𝑑θ𝑑𝑧₀ = 1/𝑘(contorno[1][2]) * (w * ρ(contorno[1][2]) * ΔH - j * Va)
bvp1 = BVProblem(modelo_marc!, bc!, [300.,500.], (0., 15.))
sol = solve(bvp1, GeneralMIRK4(), dt=0.01)
plot(sol)

