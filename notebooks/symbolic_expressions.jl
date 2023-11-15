### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ f4bbad10-7e14-11ee-0343-7da1e5ed0283
# ╠═╡ show_logs = false
begin
	cd(mktempdir())
	import Pkg
	Pkg.activate(".")
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl"),
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/SymbolicThreeBodyDecays.jl"),
		Pkg.PackageSpec("Plots"), Pkg.PackageSpec("Parameters"), Pkg.PackageSpec("JSON")
	])
	# 
	using JSON
	using Parameters
	using LinearAlgebra
	using ThreeBodyDecay
	using SymbolicThreeBodyDecays
	using SymbolicThreeBodyDecays.SymPy
	# 
	using Plots
	import Plots.PlotMeasures: mm
end

# ╔═╡ a5d25226-f465-428f-931b-64825963c255
md"""
# Spin-averated $\Lambda_c^{**+} \to \Lambda_c^+ \pi^+ \pi^-$ distribution
"""

# ╔═╡ b4c00f7e-ce88-4516-a43a-d63e534bad6e
md"""

The notebook computes the **spin-averated** angular functions for the three body decays of the exited $\Lambda_c^{**+}$ states under different hypotheses of its quantum numbers.

The spin factors are evaluated using conventions of the three-body decays introduced in the Dalitz plot decomposition paper.

The implementation in Julia comes from the [ThreeBodyDecay.jl](https://github.com/mmikhasenko/ThreeBodyDecay.jl) package.
A dispatch on SymPy object is used to provide result in a nice symbolic form. The dispatch is defined in [SymbolicThreeBodyDecay.jl](https://github.com/mmikhasenko/SymbolicThreeBodyDecays.jl) package.

The spin sums are evaluated symbolically. The resulting decay intensity distribution is a bilinear form on isobars lineshapes.

$I = |R^{1}|^2\, T_{1,1}(\theta,\zeta) + |R^{2}|^2\,T_{2,2}(\theta,\zeta) + 2 \text{Re} (R^{1*}R^{2})\,T_{1,2}(\theta,\zeta)$

where $I$ is a differential decay rate (density of the events on the Dalitz plot),
$R^{i}$ is the parametrization of the lineshape of an intermediate resonance $i$,
and $\theta$, $\zeta$ are kimanatic angles that are unambiguesly expressed using the Mandelstam invariants.
 
* The LS couplings are used for parametrizing the decay vertices. The minimal orbital angular momentum is chosen is several value are possible.
* The matrix $T$ gives only the angular factors, the break-up momentum for $p^L q^l$ is attributed to the lineshape function $R$, with $L$ being the orbital angular momentum in the decay of $\Lambda_c^{**+}$, and $l$ for the decay of an intermediate resonance.
* We incorporate **six** subchannel resonances: $f_0$, and $\rho$ in $\pi\pi$, as well as $\Sigma_c(2455)$ and $\Sigma_c(2520)$ for both $\Lambda_c \pi$ channels
* Value of $l$ is unambigues in all channels, $l = 0$ for $f_0 \to \pi\pi$, and $l=1$ for all other resonances.
* The value of $L$ in the decay of $\Lambda_c^{**+}$ is indicated in the lineshape function:

$\large R^{\text{resonance}}_{L,S}$


"""

# ╔═╡ 148e794d-3d4b-483b-8daf-b2526198f13a
md"""
## How to related $R^{\Sigma(2455)++}$ and $R^{\Sigma(2455)0}$?

The isospin Clebsch in the $\Lambda_c^{**+}$ decay is the same

```math
\begin{align}
\langle \Sigma(2245)^{++} ; \pi^- |\Lambda_c^{**+} \rangle &= \langle 1,1 ; 1,-1 |0,0 \rangle = \frac{1}{\sqrt{3}}\,,\\
\langle \Sigma(2245)^{0} ; \pi^+ |\Lambda_c^{**+} \rangle &= \langle 1,-1 ; 1,1 |0,0 \rangle = \frac{1}{\sqrt{3}}\,,\\
\end{align}
```

There is a difference treatment of the $\Sigma_c$ decays.

For $\Sigma_c^{++} (12)$ decaying to $\Lambda_c^{+}(1) \pi^+(2)$ in P-wave, the helicity coupling is related to the LS coupling with Clebshes which are the same

```math
\begin{align}
H^{\Sigma_c^{++}\to \Lambda_c^{+} \pi^+}_{\nu;\nu,0} &= 
\langle1,0;1/2,\nu|1/2,\nu\rangle\,
\langle1/2,\nu; 0,0|1/2,\nu\rangle\,
H^{\Sigma_c^{++}\to \Lambda_c^{+} \pi^+}_{P-\text{wave}}\,,\\
H^{\Sigma_c^{0}\to \pi^-\Lambda_c^{+}}_{\nu;0,\nu} &= 
\langle1,0;1/2,\nu|1/2,\nu\rangle\,
\langle0,0;1/2,\nu|1/2,\nu\rangle\,
H^{\Sigma_c^{0}\to \pi^-,\Lambda_c^{+}}_{P-\text{wave}}
\end{align}
```

The symmetry of the LS coupling under permultation is the simplest,

$H^{\Sigma_c\to \Lambda_c \pi}_{L-\text{wave}} = (-1)^L H^{\Sigma_c\to \pi \Lambda_c}_{L-\text{wave}}$

The relation gives the way to relate the lineshape functions of $\Sigma_c^{0}$ and $\Sigma_c^{++}$:

$R^{\Sigma_c^{0}} = -R^{\Sigma_c^{++}}\,.$

It holds for both $\Sigma_c(2455)$ and $\Sigma_c(2520)$.
"""

# ╔═╡ 002b5c97-f66d-4ebd-b040-b3e341675575
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto, :auto), ylim=(:auto, :auto),
    lw=1.2, lab="", colorbar=false,
    bottom_margin=3mm, left_margin=3mm)

# ╔═╡ d9c8a2e6-869e-4dbc-8c0d-8e88d69e1a4f
md"""
## Integration with Symbolic Computation
"""

# ╔═╡ cf4a76bd-a265-47be-a26c-71ebdc541838
begin
    @syms m1 m2 m3 m0
    @syms σ1 σ2 σ3
end;

# ╔═╡ dd425a12-1ffc-4eb3-9a05-64bdbfce4709
const ms = ThreeBodyMasses(m1, m2, m3; m0)

# ╔═╡ d7b53398-05a8-4bb6-92d3-bfef3991ed2b
const σs = Invariants(ms; σ1, σ2)

# ╔═╡ 2fe649e9-1f89-43ba-90f6-3243e15d85b0
function spinparity(p)
    pt = (p[2]..., p[1])
    jpv = str2jp.(pt)
    getproperty.(jpv, :j) .|> x2 |> ThreeBodyDecay.SpinTuple,
    getproperty.(jpv, :p) |> ThreeBodyDecay.ParityTuple
end

# ╔═╡ aaede4c1-1c32-4baf-b810-6884dcbe755e
begin
	Base.@kwdef struct Setup{T}
		two_js
		decay_chains::Vector{T}
	end
	function Setup(chains::Vector{Tuple{String, String, Int64}}, ifstate)
	    js, Ps = ifstate |> spinparity
	    tbs = ThreeBodySystem(ms, js)
		# 
		dcv = map(chains) do chain
			nameR, jpR, k = chain
			Rjp = str2jp(jpR)
			# 
			d0 = DecayChainLS(k, identity; two_s=Rjp.j |> x2, parity=Rjp.p, Ps, tbs)
			two_L, two_S = d0.HRk.two_ls
			L = div(two_L, 2)
			symR = sympy.Symbol(nameR * "_{$(L), $(two_S)/2}", real=true)
			d = DecayChainLS(k, σ->symR; two_s=Rjp.j |> x2, parity=Rjp.p, Ps, tbs)
		end
		return Setup(tbs.two_js, dcv)
	end
end

# ╔═╡ 3174e990-f4a2-456f-b6e6-5b810eb6c32b
function unpolarized_intensity(setup; refζs=(1, 1, 1, 1))
	# 
	@unpack two_js, decay_chains = setup
	full_amplitude = sum(itr(two_js)) do two_λs
        A = sum(decay_chains) do dc 
			amplitude(dc, σs |> StickySymTuple, two_λs .|> Sym; refζs)
		end
		abs2(A)
    end
    full_amplitude
end

# ╔═╡ 122832c5-17ab-449c-aacf-30a519f4b983
md"""
## Nice printing
"""

# ╔═╡ 8488e3e1-06e1-4800-9cc0-55a21e7deeb9
latexalign(H::Matrix) = 
"""
```math
\\small
\\begin{align}
""" *
	prod(map(eachindex(Symmetric(H))) do i
		i[1]!=i[2] ? "" :
		"T_{$(i[1]),$(i[2])} &= " * sympy.latex(H[i[1],i[2]]) * "\\,,\\\\"
	end) *	
	prod(map(eachindex(Symmetric(H))) do i
		i[1]<i[2] ?
		"T_{$(i[1]),$(i[2])} &= " * sympy.latex(H[i[1],i[2]]) * "\\,,\\\\" : ""
	end) *
	prod(map(eachindex(Symmetric(H))) do i
		i[1]>i[2] ? "T_{$(i[1]),$(i[2])} &= T_{$(i[2]),$(i[1])}" * "\\,,\\\\" : ""
	end) *
"""
\\end{align}
```
"""

# ╔═╡ a86d6bb7-7108-48b7-8597-4a9bae85b351
latexintensity(V::Vector, Isubscript="") = 
"""
```math
\\small
\\begin{equation}
{\\Large I_{$(Isubscript)}}  = 
\\begin{pmatrix}
""" *
	sympy.latex(V[1]) *
	prod(map(V[2:end]) do v
		" \\\\ " * sympy.latex(v)
	end) *	
"""
\\end{pmatrix}^{\\Large\\dagger}
\\,\\, {\\Huge [T]} \\,\\,
\\begin{pmatrix}
""" *
	sympy.latex(V[1]) *
	prod(map(V[2:end]) do v
		" \\\\ " * sympy.latex(v)
	end) *	
"""
\\end{pmatrix}
\\end{equation}
```
"""

# ╔═╡ 434a94a8-f80a-438a-91c9-f7537f386c12
printalign(H::Matrix) = Markdown.parse(latexalign(H))

# ╔═╡ 36fda437-1c28-4c57-9f44-0661e07f9aff
md"""
## Application
"""

# ╔═╡ 29787914-e344-41b9-8ce3-d72f4814dde3
final_state = ("1/2+", "0-", "0-")

# ╔═╡ d280c83e-be3f-4460-8034-1ffb8accc204
reaction = "3/2-" => final_state

# ╔═╡ dedf9cc4-68fd-4ab8-b5db-b249a8bde20e
chains = [
	("R^{\\rho(770)0}", "1-", 1), #pi+pi-
	("R^{f_0}", "0+", 1), #pi+pi-
	# 
	("R^{\\Sigma(2455)0}", "1/2+", 2), #pi-Lc+
	("R^{\\Sigma(2520)0}", "3/2+", 2), #pi-Lc+
	# 
	("R^{\\Sigma(2455)++}", "1/2+", 3), #Lc+pi+
	("R^{\\Sigma(2520)++}", "3/2+", 3), #Lc+pi+
]

# ╔═╡ 97dcb79a-5e4b-4678-9476-7ac7c11b95c8
md"""
## 1/2- decay
"""

# ╔═╡ c69a0e3c-3900-4ef1-963c-33d2d08fa1aa
md"""
## 1P states: 1/2-, 3/2-
"""

# ╔═╡ da49b68a-d349-42d7-a4ce-921fedf9a007
s1h⁻ = Setup(chains, "1/2-" => final_state)

# ╔═╡ 77faf88e-6a84-4888-9c85-268da789fae1
s3h⁻ = Setup(chains, "3/2-" => final_state)

# ╔═╡ 9657fbe0-7369-4d6e-8a10-0cc4e776f366
md"""
## 1D states: 3/2+ 5/2+
"""

# ╔═╡ 803cedb0-5de6-48a4-9f87-b4d0b91c2871
s3h⁺ = Setup(chains, "3/2+" => final_state)

# ╔═╡ 552d49c7-c05d-4106-932f-1c5ab7af140e
s5h⁺ = Setup(chains, "5/2+" => final_state)

# ╔═╡ e344aac4-ab78-4052-bb73-6ab008e3cdb5
md"""
## Save to file
"""

# ╔═╡ e0640daa-fa00-421c-99ec-be93347cb751


# ╔═╡ 19cee35a-e13c-4d9f-a658-b2436b2bea89
function wrap2dict(TV)
	T,V = TV
	d_sympy = Dict()
	d_cc = Dict()
	for i in 1:size(T,1), j in 1:size(T,2)
		d_sympy["T_$i,$j"] = sympy.srepr(T[i,j])
		d_cc["T_$i,$j"] = sympy.ccode(T[i,j])
	end
	Dict(
		"lineshapes" => sympy.srepr.(V),
		"intensity_latex" => latexintensity(V),
		"spin_averaged_crossings" => Dict(
			"latex" => latexalign(T),
			"sympy" => d_sympy,
			"ccode" => d_cc))
end

# ╔═╡ b7103bf9-f46d-4ded-b44b-3c220064324f
function writejson(path, obj)
    open(path, "w") do io
        JSON.print(io, obj, 4)
    end
end

# ╔═╡ ce51ae16-da09-4c38-a0d6-902c28725687
md"""
## Simplifications
"""

# ╔═╡ bc942632-e288-48c8-afdd-65831121954e
sign_replacements = let
	all_rotations = [(i,j,k) for (i,j,k)
		in Iterators.product(1:3,1:3,0:3) if i!=j]
	negative_rotations = filter(all_rotations) do (i,j,k)
		_wr = wr(i,j,k)
		!ispositive(_wr)
	end
	map(negative_rotations) do (i,j,k)
		symζ(wr(i,j,k)) => -symζ(wr(j,i,k))
	end
end

# ╔═╡ 0f0aa90f-c577-4025-b112-023ae5e2048b
uncombine_replacement = [
	symζ(wr(2,3,1)) => symζ(wr(2,1,1))+symζ(wr(1,3,1)),
	symζ(wr(3,1,2)) => symζ(wr(3,2,2))+symζ(wr(2,1,2)),
	symζ(wr(1,2,3)) => symζ(wr(1,3,3))+symζ(wr(3,2,3))]

# ╔═╡ 1c81d29e-246f-48d6-b3da-9d1ce4a602ca
combine_replacements = [
	symζ(wr(1,3,1)) => symζ(wr(2,3,1)) - symζ(wr(2,1,1)),
	symζ(wr(2,1,2)) => symζ(wr(3,1,2)) - symζ(wr(3,2,2)),
	symζ(wr(3,2,3)) => symζ(wr(1,2,3)) - symζ(wr(1,3,3))]

# ╔═╡ 9f9a1ddc-c8ec-47d3-84a3-c4c3c4e2efbf
function intensity_matrix(setup)
	I_symm = unpolarized_intensity(setup; refζs=(1, 1, 1, 1));
	lineshapes = map(d->d.Xlineshape(nothing), setup.decay_chains)
	H = sympy.hessian(I_symm, lineshapes)
	H_fully_simp = map(eachindex(Symmetric(H))) do i
		i[1]>i[2] ? Sym(0) : 
		simplify(H[i].doit())
	end |> Symmetric;
	H_normalized = H_fully_simp.
			xreplace(sign_replacements |> Dict).
			xreplace(combine_replacements |> Dict)#.
			# xreplace(unsafe_uncombine_replacement[2:3] |> Dict);
	return H_normalized, lineshapes
end

# ╔═╡ 26d07dd7-ea69-41f1-a959-940581497f72
@time T1h⁻ = intensity_matrix(s1h⁻);

# ╔═╡ 64f27f9c-5555-4408-a650-f36559e3fe2b
latexintensity(T1h⁻[2], "1/2-") |> Markdown.parse

# ╔═╡ be83f923-a9bc-4aa1-ac66-21ebd69eedc7
printalign(T1h⁻[1])

# ╔═╡ c1db0f6c-d47d-474b-97fa-2582a02c35a8
@time T3h⁻ = intensity_matrix(s3h⁻);

# ╔═╡ 4d14d9b2-dd03-4697-af66-a2667c0b4509
latexintensity(T3h⁻[2], "3/2-") |> Markdown.parse

# ╔═╡ ff1cfcff-661b-42d3-97d3-0d54bc7552a5
printalign(T3h⁻[1])

# ╔═╡ 31aa74d1-a5c5-44fe-905c-7e75da3e8195
@time T3h⁺ = intensity_matrix(s3h⁺);

# ╔═╡ 11cb20ea-e9e9-4914-ae97-e484396008c8
latexintensity(T3h⁺[2], "3/2+") |> Markdown.parse

# ╔═╡ e992cb27-c7e2-416d-b02e-77e3b3eb2aed
printalign(T3h⁺[1])

# ╔═╡ 4a2e4f68-41f7-4273-a30d-26bd052f0ab4
@time T5h⁺ = intensity_matrix(s5h⁺);

# ╔═╡ 4390fe7b-8fb2-4793-9743-6717eacc1e9d
latexintensity(T5h⁺[2], "5/2+") |> Markdown.parse

# ╔═╡ 3e9c2854-2e4b-4bdf-8dba-046a277ab078
printalign(T5h⁺[1])

# ╔═╡ ee898d12-e33c-49c4-a3bd-e5ba91be96b3
computation_result = Dict("1/2-"=>wrap2dict(T1h⁻),
	 "3/2-"=>wrap2dict(T3h⁻),
	 "3/2+"=>wrap2dict(T3h⁺),
	 "5/2+"=>wrap2dict(T5h⁺));

# ╔═╡ 747f4e33-44a9-48fb-8297-fe6f83a49518
writejson(joinpath(@__DIR__, "LcXX2Lcpipi.json"), computation_result)

# ╔═╡ 8151b9a8-2d75-4672-b290-2cdeb0c89cb0
unsafe_uncombine_replacement = [
	symζ(wr(2,3,0)) => 2PI-symζ(wr(3,1,0))-symζ(wr(1,2,0)),
	symζ(wr(3,1,0)) => 2PI-symζ(wr(1,2,0))-symζ(wr(2,3,0)),
	symζ(wr(1,2,0)) => 2PI-symζ(wr(2,3,0))-symζ(wr(3,1,0))
]

# ╔═╡ Cell order:
# ╟─a5d25226-f465-428f-931b-64825963c255
# ╟─b4c00f7e-ce88-4516-a43a-d63e534bad6e
# ╟─148e794d-3d4b-483b-8daf-b2526198f13a
# ╠═f4bbad10-7e14-11ee-0343-7da1e5ed0283
# ╠═002b5c97-f66d-4ebd-b040-b3e341675575
# ╟─d9c8a2e6-869e-4dbc-8c0d-8e88d69e1a4f
# ╠═cf4a76bd-a265-47be-a26c-71ebdc541838
# ╠═dd425a12-1ffc-4eb3-9a05-64bdbfce4709
# ╠═d7b53398-05a8-4bb6-92d3-bfef3991ed2b
# ╠═2fe649e9-1f89-43ba-90f6-3243e15d85b0
# ╠═aaede4c1-1c32-4baf-b810-6884dcbe755e
# ╠═3174e990-f4a2-456f-b6e6-5b810eb6c32b
# ╠═9f9a1ddc-c8ec-47d3-84a3-c4c3c4e2efbf
# ╟─122832c5-17ab-449c-aacf-30a519f4b983
# ╠═8488e3e1-06e1-4800-9cc0-55a21e7deeb9
# ╠═a86d6bb7-7108-48b7-8597-4a9bae85b351
# ╠═434a94a8-f80a-438a-91c9-f7537f386c12
# ╟─36fda437-1c28-4c57-9f44-0661e07f9aff
# ╠═29787914-e344-41b9-8ce3-d72f4814dde3
# ╠═d280c83e-be3f-4460-8034-1ffb8accc204
# ╠═dedf9cc4-68fd-4ab8-b5db-b249a8bde20e
# ╟─97dcb79a-5e4b-4678-9476-7ac7c11b95c8
# ╟─c69a0e3c-3900-4ef1-963c-33d2d08fa1aa
# ╠═da49b68a-d349-42d7-a4ce-921fedf9a007
# ╠═26d07dd7-ea69-41f1-a959-940581497f72
# ╠═64f27f9c-5555-4408-a650-f36559e3fe2b
# ╠═be83f923-a9bc-4aa1-ac66-21ebd69eedc7
# ╠═77faf88e-6a84-4888-9c85-268da789fae1
# ╠═c1db0f6c-d47d-474b-97fa-2582a02c35a8
# ╠═4d14d9b2-dd03-4697-af66-a2667c0b4509
# ╠═ff1cfcff-661b-42d3-97d3-0d54bc7552a5
# ╟─9657fbe0-7369-4d6e-8a10-0cc4e776f366
# ╠═803cedb0-5de6-48a4-9f87-b4d0b91c2871
# ╠═31aa74d1-a5c5-44fe-905c-7e75da3e8195
# ╠═11cb20ea-e9e9-4914-ae97-e484396008c8
# ╠═e992cb27-c7e2-416d-b02e-77e3b3eb2aed
# ╠═552d49c7-c05d-4106-932f-1c5ab7af140e
# ╠═4a2e4f68-41f7-4273-a30d-26bd052f0ab4
# ╠═4390fe7b-8fb2-4793-9743-6717eacc1e9d
# ╠═3e9c2854-2e4b-4bdf-8dba-046a277ab078
# ╟─e344aac4-ab78-4052-bb73-6ab008e3cdb5
# ╠═e0640daa-fa00-421c-99ec-be93347cb751
# ╠═19cee35a-e13c-4d9f-a658-b2436b2bea89
# ╠═b7103bf9-f46d-4ded-b44b-3c220064324f
# ╠═ee898d12-e33c-49c4-a3bd-e5ba91be96b3
# ╠═747f4e33-44a9-48fb-8297-fe6f83a49518
# ╟─ce51ae16-da09-4c38-a0d6-902c28725687
# ╠═bc942632-e288-48c8-afdd-65831121954e
# ╠═0f0aa90f-c577-4025-b112-023ae5e2048b
# ╠═1c81d29e-246f-48d6-b3da-9d1ce4a602ca
# ╠═8151b9a8-2d75-4672-b290-2cdeb0c89cb0
