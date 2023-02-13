### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 1bd1acfa-9c9e-11ed-3779-95709ce9906b
begin
	import Pkg
	Pkg.activate("/home/reginaldo/Insync/Trabalho/IFMG/IFMG_ARCOS/codigosemdesenvolvimento/EDOMARCELO/edo_marcelo_jl")
	Pkg.add("ModelingToolkit")
end

# ╔═╡ c0fca042-8f04-4bd7-abc6-178f247129e6
using Plots, ModelingToolkit, DifferentialEquations, Symbolics

# ╔═╡ fcd6d4b9-eb65-494a-ae69-f54d49d285bb
md"""
# Definindo parâmetros dependentes da temperatura $\theta$ (°C)

## Coeficiente de condutividade térmica $k$.

$$k(\theta) =  \left\{ 
\begin{array}{ll}
54. -3.33 T^2 & \textrm{; $20°C \leq \theta \leq 800°C$}\\
27. & \textrm{; $800°C < \theta \leq 1515°C$} \end{array} 
\right.$$
"""

# ╔═╡ cb8e5509-d4ad-4b45-a9e6-101471ddfe29
k(θ) = if θ >= 20. && θ <= 800.
	54 - 3.33e-2 * θ
elseif θ > 800. && θ <= 1515.
	27.3
end

# ╔═╡ 26174a5c-11be-40d6-8896-0aa30fd0e270
md"""
## Capacidade térmica volumétrica $c_v$.

$$c_v(\theta) =  \left\{ 
\begin{array}{ll}
425.0 + 7.73 . 10^{-1} T - 1.69 . 10^{-3} T^2 + 2.22 . 10^{-6} T^3 & \textrm{; $20°C \leq \theta \leq 600°C$}\\
666. + \frac{13002}{738 - T} & \textrm{; $600°C < \theta \leq 735°C$}\\
545. + \frac{17820}{T - 731} & \textrm{; $735°C < \theta \leq 900°C$}\\
650. & \textrm{; $900°C < \theta \leq 1515°C$}
\end{array} 
\right.$$
"""

# ╔═╡ 8235d805-4703-4ca2-8af7-c2f3d5ac770b
cv(θ) = if θ >= 20. && θ <= 600.
	425. + 7.73e-1 * θ - 1.69e-3 * θ^2 + 2.22e-6 * θ^3
elseif θ > 600. && θ <= 735.
	666. + 13002. / (738. - θ)
elseif θ > 735. && θ <= 900.
	525. + 17820. / (θ - 731.)
elseif θ > 900. && θ <= 1515.
	650.
end

# ╔═╡ 3bc6a4d6-cee7-4c5c-81af-de2337615dcf
md"""
## Densidade $\rho$.

$$\rho(\theta) =  \left\{ 
\begin{array}{ll}
-0.3373 T + 7871. & \textrm{; $25°C \leq \theta \leq 689°C$}\\
0.1226 T + 7556. & \textrm{; $689°C < \theta \leq 853°C$}\\
-0.5575 T + 8158. & \textrm{; $853°C < \theta \leq 1515°C$}
\end{array} 
\right.$$
"""

# ╔═╡ baa33975-0c09-4092-b15d-e668f7531dae
ρ(θ) = if θ >= 25. && θ <= 689.
	-0.3373 * θ + 7871.
elseif θ > 689. && θ <= 853.
	0.1226 * θ + 7556.
elseif θ > 853. && θ <= 1515.
	-0.5575 * θ + 8158.
end

# ╔═╡ 8652eeea-e54c-4c08-9dd4-1ab8c9931137
md"""
## Coeficiente de troca de calor $h$.

$$\rho(\theta) =  \left\{ 
\begin{array}{ll}
1175. & \textrm{; $0 mm \leq h \leq 9 mm$}\\
115.  & \textrm{; $9 mm < h \leq 15 mm$}
\end{array} 
\right.$$
"""

# ╔═╡ dee5c6f6-b527-4c5c-97c3-c785d5ca6aba
h(z) = if z >= 0. && z <= 9.
	1175.
elseif z > 9. && z <= 15.
	115.
end

# ╔═╡ f6ff5870-9720-4f4f-851c-b8567a2825db
begin
	@register_symbolic k(θ)
	@register_symbolic cv(θ)
	@register_symbolic ρ(θ)
	@register_symbolic h(z)
end

# ╔═╡ c4c27d96-d28b-49fe-bd72-2ceb885fe553
k(Symbolics.Num(200.))

# ╔═╡ b1289e5f-7939-410a-865b-7e49e2623211
md"""
# Modelo matemático

## Distribuição de temperatura em eletrodo ao longo de seu comprimento energizado na soldagem subaquática com arame tubular.

$$\frac{d}{dz}\left( k(\theta) \frac{d\theta}{dz}\right) + \omega \rho(\theta) 
\frac{d(c_v(\theta) . \theta)}{dz} - h(z).(\theta - \theta_\infty) = -j^2 r$$
"""

# ╔═╡ ef4684d2-87e7-4227-a3cc-fff4189c9946
@variables z θ(..);

# ╔═╡ bf28a50d-ec9d-4f74-b99f-a31c6b4544d4
@constants ω = 0.075 j = 99.47 va = 20. ΔH = 120. r = 5.62e-6 θ∞ = 25.;

# ╔═╡ 87cd9450-950e-462b-a2f3-604b2bd8f1c0
dz = Differential(z);

# ╔═╡ da254df5-5819-400e-b41f-71170fe9d6d8
d2z = Differential(z)^2;

# ╔═╡ b77c6af9-f89c-46b7-a888-4085766e097f
parc1 = k(θ(z)) * d2z(θ(z))

# ╔═╡ d6de7e4f-be2a-4cff-8e0b-2627e9e0a20e
parc2 = ω * ρ(θ(z)) * (θ(z) + dz(θ(z)))

# ╔═╡ 85234c73-246e-4bb1-8621-835ff3d3a1f7
parc3 = -h(z) * (θ(z) - θ∞)

# ╔═╡ 2ec8b81e-0e8e-4d55-8452-50985a0895e6
parc4 = -j^2 * r

# ╔═╡ d5c7d169-fd82-4572-b3a4-2c9cf1499f52
eq = [parc1 + parc2 + parc3 ~ parc4]

# ╔═╡ 1271e7ab-3d4b-416a-8916-45a3f3496183
cdct = [
	θ(0.) ~ 1515.,
	θ(15.) ~ 158.,
	dz(θ(z)) ~ 1/k(1515.)*(ω * ρ(1515.) * ΔH - j * va)
]

# ╔═╡ 4f2e444c-3a32-46c3-9fa7-b31aa00779c9
dominio = [
	θ ∈ (158., 1515.),
	z ∈ (0., 9.)
]

# ╔═╡ b39d8a93-461e-41bd-95dc-5dd731f1f94e
#@named pde_system = PDESystem(eq,cdct,dominio,[z],[θ])
@named ode_system = ODESystem(eq)

# ╔═╡ 804b6f50-da0e-4125-8146-f68b2acc7572
sol = symbolic_discretize(pde_system)

# ╔═╡ Cell order:
# ╠═1bd1acfa-9c9e-11ed-3779-95709ce9906b
# ╠═c0fca042-8f04-4bd7-abc6-178f247129e6
# ╟─fcd6d4b9-eb65-494a-ae69-f54d49d285bb
# ╠═cb8e5509-d4ad-4b45-a9e6-101471ddfe29
# ╟─26174a5c-11be-40d6-8896-0aa30fd0e270
# ╠═8235d805-4703-4ca2-8af7-c2f3d5ac770b
# ╟─3bc6a4d6-cee7-4c5c-81af-de2337615dcf
# ╠═baa33975-0c09-4092-b15d-e668f7531dae
# ╟─8652eeea-e54c-4c08-9dd4-1ab8c9931137
# ╠═dee5c6f6-b527-4c5c-97c3-c785d5ca6aba
# ╠═f6ff5870-9720-4f4f-851c-b8567a2825db
# ╠═c4c27d96-d28b-49fe-bd72-2ceb885fe553
# ╠═b1289e5f-7939-410a-865b-7e49e2623211
# ╠═ef4684d2-87e7-4227-a3cc-fff4189c9946
# ╠═bf28a50d-ec9d-4f74-b99f-a31c6b4544d4
# ╠═87cd9450-950e-462b-a2f3-604b2bd8f1c0
# ╠═da254df5-5819-400e-b41f-71170fe9d6d8
# ╠═b77c6af9-f89c-46b7-a888-4085766e097f
# ╠═d6de7e4f-be2a-4cff-8e0b-2627e9e0a20e
# ╠═85234c73-246e-4bb1-8621-835ff3d3a1f7
# ╠═2ec8b81e-0e8e-4d55-8452-50985a0895e6
# ╠═d5c7d169-fd82-4572-b3a4-2c9cf1499f52
# ╠═1271e7ab-3d4b-416a-8916-45a3f3496183
# ╠═4f2e444c-3a32-46c3-9fa7-b31aa00779c9
# ╠═b39d8a93-461e-41bd-95dc-5dd731f1f94e
# ╠═804b6f50-da0e-4125-8146-f68b2acc7572
