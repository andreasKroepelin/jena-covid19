### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 847ba008-11ae-11eb-1055-e39d77d59e4c
using HTTP, CSV, DataFrames, TypedTables, Statistics, GLM, Plots, PlutoUI, PrettyTables, Dates

# ╔═╡ ec4f0f72-19ee-11eb-1d75-e9a73bbdd49c
md"""
# Development of the Sars-CoV-2 pandemic in the city of Jena
"""

# ╔═╡ dcca0fc6-5090-11eb-3257-3dfa9f1178de
plotly()

# ╔═╡ 0979085c-19ef-11eb-38b8-032ae5772061
md"## Retrieve raw data"

# ╔═╡ c6da68b2-11ae-11eb-2bdc-3d7da6237852
csv = CSV.File(HTTP.get("https://opendata.jena.de/dataset/2cc7773d-beba-43ad-9808-a420a67ffcb3/resource/d3ba07b6-fb19-451b-b902-5b18d8e8cbad/download/corona_erkrankungen_jena.csv").body,
	normalizenames = true
	);

# ╔═╡ e9eaed36-11ae-11eb-21a9-2b5d76dd2d60
raw = dropmissing!(DataFrame(csv), disallowmissing=true);

# ╔═╡ ffd22e44-1707-11eb-3f65-b16a0dacc6f0
raw.zeit_ = unix2datetime.(raw.zeit);

# ╔═╡ bd02c648-19fd-11eb-3e30-f98ecf4304fe
md"""Last data update: **$(Dates.format(last(raw.zeit_), "Y u d, H:mm"))**"""

# ╔═╡ 0e160200-19fd-11eb-1848-215741cdc13e
names(raw)

# ╔═╡ 740fb120-19ef-11eb-357c-1dd462cdf1ba
md"## Restructure data into equidistant time steps"

# ╔═╡ 3db6654c-1709-11eb-2433-b1a4e12530d6
data = let
	time = let
		stop = last(raw.zeit_)
		start = Dates.tonext(floor(first(raw.zeit_), Day); step=Minute(1), same=true) do d
			hour(d) == hour(stop) &&
			minute(d) == minute(stop)
		end 
		start : Day(1) : stop
	end
	
	corresp_idcs = [argmin(abs.(raw.zeit_ .- t)) for t in time]
	
	Table(
		time = time,
		cases = [raw.erkrankte[i] for i in corresp_idcs],
		active = [raw.aktive_faelle[i] for i in corresp_idcs],
		recovered = [raw.genesene[i] for i in corresp_idcs],
		new_cases = [raw.neu_erkrankte[i] for i in corresp_idcs],
		dead = [raw.tote[i] for i in corresp_idcs]
	)
end;

# ╔═╡ f93376d2-7d2c-11eb-33b6-53b201381a37
HTML() do io
	print(io, """<div style="height: 300px; overflow: auto;"> """)
	pretty_table(io, data, backend=:html, standalone=false)
	print(io, """</div> """)
end

# ╔═╡ 45eca0fc-19f2-11eb-18a1-4f9d0791c900
md"### Total number of cases"

# ╔═╡ 0b779326-170b-11eb-1602-7dbe1e277b45
plot(data.time, data.cases, label=:none)

# ╔═╡ 6885f3e8-19f2-11eb-0b89-4f47857177f1
md"### Cases and recovered"

# ╔═╡ 94fa8d5c-19f3-11eb-025f-51738530f25d
cases_recovered_lag = let
	function diff(l)
		a = @view data.cases[1:end-l]
		b = @view data.recovered[l+1:end]
		(a .- b) .^ 2 |> mean
	end
	argmin(diff.(0:80))
end;

# ╔═╡ 9d148f70-19f2-11eb-2c10-2104da42ddbd
md"""
Lead recovered curve by: $(@bind recov_lead Slider(0:50, show_value=true)) days

By minimising the squared difference between the curves, I find that the recovered curve lags behind the cases curve by **$cases_recovered_lag days**.
This implies a recovery time of **$cases_recovered_lag days**.
"""

# ╔═╡ c2df3446-3bce-11eb-1517-697d0f346f6f
md"### Dead"

# ╔═╡ d457eedc-3bce-11eb-18c1-337ec36123c7
plot(data.time, data.dead, label=:none)

# ╔═╡ 97b3fe34-19fd-11eb-1d31-cb439f8ed8e0
md"### New cases per day"

# ╔═╡ 9f2cd460-19fd-11eb-3de8-457f73cd390f
bar(data.time, data.new_cases, label=:none)

# ╔═╡ d943da52-19f5-11eb-2472-7356550d3081
md"### Active cases"

# ╔═╡ f02739e2-1768-11eb-176d-416bde3d4b1a
plot(data.time, data.active, label=:none)

# ╔═╡ af8c6318-19f6-11eb-2d3b-31b5b159a98a
md"""
Show incidence over $(@bind n Slider(1:100, default=7, show_value=true)) days.

Normalise number onto period of $(@bind norm_incidence_n Slider(1:100, default=7, show_value=true)) days? $(@bind norm_incidence CheckBox())
"""

# ╔═╡ 7c155ae2-19f6-11eb-1d19-d96b7ab4e155
md"### $n days incidence per 100 000 residents"

# ╔═╡ 1fa8a168-170c-11eb-2098-2f34033078ef
let
	population = 108_127
	last_n_days = [Float64(data.cases[i] - data.cases[i-n]) for i in n+1:size(data)[1]]
	if norm_incidence
		last_n_days ./= n / norm_incidence_n
	end
	plot(data.time[n+1:end], last_n_days ./ population .* 100_000, label=:none)
end

# ╔═╡ eef68ee6-19f8-11eb-1ccd-d300689ecd65
md"### Fit exponential model to ..."

# ╔═╡ 6cc7d9fa-7d2b-11eb-05ab-7b22ed4a86a4
modelled_series = :cases

# ╔═╡ 89240fe8-19f9-11eb-2324-a7a1fb9b6398
md"""
Exponential model between $(@bind t0 DateField(default=Dates.now()-Day(7))) and $(@bind t1 DateField(default=Dates.now()))
"""

# ╔═╡ 06ed628e-1a05-11eb-2dde-8b0414c4cce1
struct ExpModel
	values::Vector{Union{Missing, Float64}}
	double_time::Float64
	R²::Float64
end

# ╔═╡ 7d347f06-170e-11eb-3618-cfd9dbc0eeb8
function exp_model(series::Symbol, t0, t1)
	days_since_t0(orig_time) = (orig_time - t0) / Millisecond(Day(1))
	time_selection = t0 .<= data.time .<= t1
	fit_tab = Table(
		values = getproperty(data, series)[time_selection],
		time = days_since_t0.(data.time[time_selection])
	)
	
	model = lm(@formula(log(values) ~ 1 + time), fit_tab)
	
	predict_tab = Table(time = days_since_t0.(data.time))
	values = predict(model, predict_tab)
	values[(data.time .< t0) .| (data.time .> t1)] .= missing
	
	ExpModel(
		exp.(values),
		log(2) / coef(model)[2],
		r2(model)
	)
end

# ╔═╡ fc557a84-19f8-11eb-3b57-03238b095c9d
exp_model_cases = exp_model(modelled_series, t0, t1);

# ╔═╡ f9c3fc86-19f9-11eb-1156-e7bdc755a3a2
md"""
This model implies that numbers **double every $(round(Int, exp_model_cases.double_time)) days** in the specified period.

Coefficient of determination ``R^2`` of that regression: $(exp_model_cases.R²)
"""

# ╔═╡ f507f098-170c-11eb-0471-d9c6fc16defa
begin
	plot(data.time, [getproperty(data, modelled_series) exp_model_cases.values], label=["real" "model"], legend=:topleft)
	plot!([t0, t1], [0, 0], lw=5, label="region of regression")
end

# ╔═╡ de29ac72-1775-11eb-2493-fbe9c4210288
lead(xs, n) = [xs[n+1:end]..., fill(NaN, n)...]

# ╔═╡ d629204a-1773-11eb-166d-871176ee1799
plot(
	data.time,
	[data.cases lead(data.recovered, recov_lead)],
	label=["cases" """recovered$(recov_lead > 0 ? " (leading by $recov_lead days)" : "")"""],
	legend=:bottomright
)

# ╔═╡ Cell order:
# ╟─ec4f0f72-19ee-11eb-1d75-e9a73bbdd49c
# ╠═847ba008-11ae-11eb-1055-e39d77d59e4c
# ╟─dcca0fc6-5090-11eb-3257-3dfa9f1178de
# ╟─0979085c-19ef-11eb-38b8-032ae5772061
# ╠═c6da68b2-11ae-11eb-2bdc-3d7da6237852
# ╠═e9eaed36-11ae-11eb-21a9-2b5d76dd2d60
# ╠═ffd22e44-1707-11eb-3f65-b16a0dacc6f0
# ╟─bd02c648-19fd-11eb-3e30-f98ecf4304fe
# ╠═0e160200-19fd-11eb-1848-215741cdc13e
# ╟─740fb120-19ef-11eb-357c-1dd462cdf1ba
# ╠═3db6654c-1709-11eb-2433-b1a4e12530d6
# ╟─f93376d2-7d2c-11eb-33b6-53b201381a37
# ╟─45eca0fc-19f2-11eb-18a1-4f9d0791c900
# ╠═0b779326-170b-11eb-1602-7dbe1e277b45
# ╟─6885f3e8-19f2-11eb-0b89-4f47857177f1
# ╟─d629204a-1773-11eb-166d-871176ee1799
# ╟─9d148f70-19f2-11eb-2c10-2104da42ddbd
# ╠═94fa8d5c-19f3-11eb-025f-51738530f25d
# ╟─c2df3446-3bce-11eb-1517-697d0f346f6f
# ╠═d457eedc-3bce-11eb-18c1-337ec36123c7
# ╟─97b3fe34-19fd-11eb-1d31-cb439f8ed8e0
# ╠═9f2cd460-19fd-11eb-3de8-457f73cd390f
# ╟─d943da52-19f5-11eb-2472-7356550d3081
# ╠═f02739e2-1768-11eb-176d-416bde3d4b1a
# ╟─7c155ae2-19f6-11eb-1d19-d96b7ab4e155
# ╠═1fa8a168-170c-11eb-2098-2f34033078ef
# ╟─af8c6318-19f6-11eb-2d3b-31b5b159a98a
# ╟─eef68ee6-19f8-11eb-1ccd-d300689ecd65
# ╠═6cc7d9fa-7d2b-11eb-05ab-7b22ed4a86a4
# ╟─89240fe8-19f9-11eb-2324-a7a1fb9b6398
# ╟─f9c3fc86-19f9-11eb-1156-e7bdc755a3a2
# ╠═fc557a84-19f8-11eb-3b57-03238b095c9d
# ╠═f507f098-170c-11eb-0471-d9c6fc16defa
# ╟─06ed628e-1a05-11eb-2dde-8b0414c4cce1
# ╠═7d347f06-170e-11eb-3618-cfd9dbc0eeb8
# ╟─de29ac72-1775-11eb-2493-fbe9c4210288
