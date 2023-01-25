### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ 1bd1acfa-9c9e-11ed-3779-95709ce9906b
begin
	import Pkg
	Pkg.activate("/home/reginaldo/Insync/Trabalho/IFMG/IFMG_ARCOS/codigosemdesenvolvimento/EDOMARCELO/edo_marcelo_jl")
   	Pkg.instantiate()
end

# ╔═╡ c0fca042-8f04-4bd7-abc6-178f247129e6
using Plots

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
𝑘(θ) = if θ >= 20 && θ <= 800
	54 - 3.33e-2 * θ
elseif θ > 800 && θ <= 1515.
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
𝒉(𝒛) = if 𝒛 >= 0. && 𝒛 <= 9.
	1175.
elseif 𝒛 > 9. && 𝒛 <= 15.
	115.
end

# ╔═╡ Cell order:
# ╠═1bd1acfa-9c9e-11ed-3779-95709ce9906b
# ╠═c0fca042-8f04-4bd7-abc6-178f247129e6
# ╟─fcd6d4b9-eb65-494a-ae69-f54d49d285bb
# ╠═cb8e5509-d4ad-4b45-a9e6-101471ddfe29
# ╟─26174a5c-11be-40d6-8896-0aa30fd0e270
# ╠═8235d805-4703-4ca2-8af7-c2f3d5ac770b
# ╠═3bc6a4d6-cee7-4c5c-81af-de2337615dcf
# ╠═baa33975-0c09-4092-b15d-e668f7531dae
# ╟─8652eeea-e54c-4c08-9dd4-1ab8c9931137
# ╠═dee5c6f6-b527-4c5c-97c3-c785d5ca6aba
