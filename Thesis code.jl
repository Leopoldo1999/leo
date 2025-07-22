using JLD2, FileIO
using CairoMakie
using Printf
using Makie
using Dierckx
using Glob
using SpecialFunctions
using LinearAlgebra
using Distributions

println("Insert the mass of the MBH in solar mass, i.e. = 1.0e4")
CMBH = parse(Float64, readline())

if CMBH < 1.0e6
    @printf("Mass of MBH: %.1f\n", CMBH)
else
    println("Mass of MBH: ", CMBH)
end

println("Which method do you want to use?")
println("1: Explore snapshot by snapshot")
println("2: Sampling from heatmap (integrate over the whole time-window)")
method = parse(Int, readline())

final_output = load("summary_2_$(CMBH)_4.8e-7_50x100.jld2")
params = final_output["params"]
aux = final_output["aux"]

time = final_output["t"]/params["Gyr"]
TDE_rate = final_output["captures"][1] * params["Gyr"] * 1e-9 #rates in 1/yr
TDE_binary = final_output["captures"][2] * params["Gyr"] * 1e-9 #rates in 1/yr

TDE_norm=TDE_rate/maximum(TDE_rate) #normalize rates
TDE_binary_norm=TDE_binary/maximum(TDE_binary) #normalize rates

fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], xlabel=L"$t$ [Gyr]", ylabel=L"\dot{N} / \dot{N}_{\max}")

lines!(ax, time, TDE_norm, label="TDEs rate")
lines!(ax, time, TDE_binary_norm, label="Tidal binary rate")

axislegend(ax, position=:rb)

display(fig)

function trapz(x, y)
    sum( 0.5 * (y[1:end-1] .+ y[2:end]) .* diff(x) )
end

t=[]
folder = "snapshots_2_$(params["M"])_4.8e-7"
basename_prefix = "snap_grid_$(length(aux["R_val"]))x$(length(aux["s_val"]))_"
files = glob(basename_prefix * "*", folder)

indices_snapshot = 1:length(files)-1 #explore how many snapshot there are

for i in indices_snapshot
    push!(t,time[i*40])
end

t_new=vcat(0,t)

println("Insert the time t1 in Gyr:")
t1 = parse(Float64, readline())
println("Insert the time t2 in Gyr:")
t2 = parse(Float64, readline())
max_ind = findall(x -> t1 <= x <= t2, t)

fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], xlabel=L"$t$ [Gyr]", ylabel=L"\dot{N}\ [1/yr] ")

lines!(ax, time, TDE_binary, label="TDEs binary")

vlines!(ax, [t1, t2], color=:red, linewidth=2, linestyle=:dash)

axislegend(ax, position=:rb)

save("Selected timespan $(params["M"])_0.1AU.pdf", fig)

E = params["E0"] * (exp.(aux["s_val"]) .- 1.0)

function besseli_leo(x,y)
    return besseli(0, x)./besseli(1, y)
end

q_t=[]
f_t=[]
dN_dEdt=[]

for i in max_ind
    snapshot_file = "snapshots_2_$(params["M"])_4.8e-7/snap_grid_$(length(aux["R_val"]))x$(length(aux["s_val"]))_$(@sprintf("%03d", i)).jld2"
    snap = load(snapshot_file)
    push!(dN_dEdt,snap["diff_rate"][2] * params["Gyr"])
    push!(q_t, snap["q"][2])
    push!(f_t, snap["f_e"][2])
end

dN_dEdt = Matrix(reduce(hcat, dN_dEdt)')
q_t = Matrix(reduce(hcat, q_t)')


N_E = length(E)
r_LC = final_output["params"]["rlc"][2]

N_rp = 100
r_p = range(0, r_LC; length=N_rp)
dx = step(r_p)

dNdrpdEdt = zeros(N_E, length(max_ind), N_rp)
dN_dEdt_2=dN_dEdt*0.0

for i in 1:N_E
    q = q_t[:,i]
    f = [ f_t[j][aux["Rlc_idx"][2][i], i] for j in 1:length(max_ind)]
    dNdrp_matrix = zeros(length(max_ind), length(r_p))
    for j in 1:length(max_ind)
        dN_dEdt_2[:,i]=4*pi^2*aux["l2c"][i]*q.*(q.^2+q.^4).^-0.25 .* f
        arg = @. 2 * sqrt.(r_p ./ (q[j] * r_LC))
        norm = 2 / sqrt(q[j])
        if q[j] < 0.01
            dN_drp = zeros(length(r_p))
            x0 = r_LC
            ix = findall(abs.(r_p .- x0) .< dx/2)
            dN_drp[ix] .= 1 / dx
        elseif q[j] > 10
            dN_drp = fill(1/r_LC, length(r_p))
        else
            dN_drp = (1/(r_LC*sqrt(q[j]))) .* @.besseli_leo(arg, norm)
        end
        dNdrp_matrix[j, :] .= dN_drp
    end
    dNdrpdEdt[i, :, :] = dN_dEdt_2[:,i] .* dNdrp_matrix
    
end

dNdrpdE = [trapz(t_new[max_ind], dNdrpdEdt[i, :, j]) for i in 1:length(E), j in 1:length(r_p)]

fig = Figure(size=(800, 600))
ax = Axis(fig[1,1], xlabel=L"Specific\ energy\ in\ internal\ units", ylabel=L"r_{p}/r_{lc}", limits=(-200,10,nothing,nothing), title="Distribution of E and rₚ between t = $t1 and t = $t2 Gyr")
hm = heatmap!(ax, -(E.-params["phis0"]), r_p/r_LC, dNdrpdE; colormap=:viridis, colorrange=(0, 0.5*maximum(dNdrpdE)))
Colorbar(fig[1,2], hm,  label = L"\frac{dN}{dE \, dr_p}")

save("Heatmap $(params["M"])_0.1AU.pdf", fig)

dt = diff(t_new*1e9)
sub_TDE_binary=TDE_binary[1:40:end]
lambas=sub_TDE_binary[1:length(files)-1].*dt
lamba2=trapz(t_new[max_ind],sub_TDE_binary[max_ind])

realizations= 1:1
TDE_final=[]
for j in realizations
    Total_TDE=[]
    for i in lambas
        y=rand(Poisson(i))
        push!(Total_TDE,y)
    end
    push!(TDE_final,Total_TDE[max_ind])
end

snap1=[]
time_snap=[]
	
for i in 1:length(max_ind)
	if TDE_final[1][i]!=0
		for j in 1:TDE_final[1][i]
			push!(snap1, max_ind[i])
			push!(time_snap, t_new[max_ind[i]])
		end
	end
end

using Random, Interpolations

GMsun_km = 1.3271244e11       # km^3/s^2,
AU_in_km = 1.495978707e8      # km to AU
N_events=sum(TDE_final[1])

if method == 1

	accepted_energies_final=[]
	for j in realizations
		accepted_energies = []
		for i in 1:length(max_ind)
			itp = interpolate((E,), dN_dEdt_2[i,:], Gridded(Linear()))
			
			p_max = maximum(dN_dEdt_2[i,:])
		
			N_TDE = TDE_final[j][i]
		
			E_min = minimum(E)
			E_max = maximum(E)
		
			Ei = Float64[]
			
			while length(Ei) < N_TDE
				E_trial = rand() * (E_max - E_min) + E_min
				y = rand() * p_max
				if y < itp(E_trial)
					push!(Ei, -(E_trial-params["phis0"]))
				end
			end
		
			push!(accepted_energies, Ei)
		end
		push!(accepted_energies_final, accepted_energies)
	end

	acepted_rp_total=[]
	acepted_B_total=[]

	for k in realizations
		accepted_rp = []
		accepted_B = []
		for i in max_ind
			snapshot_file = "snapshots_2_$(params["M"])_4.8e-7/snap_grid_$(length(aux["R_val"]))x$(length(aux["s_val"]))_$(@sprintf("%03d", i)).jld2"
			snap = load(snapshot_file)
			energies = accepted_energies_final[k][i-minimum(max_ind)+1]
			q_E = snap["q"][2]
			itp_qE = interpolate((E,), q_E, Gridded(Linear()))
			rp= Float64[]
			B= Float64[]
			for j in energies
				q=itp_qE[-(j)+params["phis0"]]
				arg = @. 2 * sqrt.(r_p ./ (q * r_LC))
				norm = 2 / sqrt(q)
				if q < 0.01
					push!(rp,r_LC)
					push!(B, 1)
				elseif q > 10
					push!(rp, rand() * r_LC)
					push!(B, r_LC / (rand() * r_LC))
				else
					dN_drp = (1/(r_LC*sqrt(q))) .* @.besseli_leo(arg, norm)
					while true  
						itp_rp = interpolate((collect(r_p),), dN_drp, Gridded(Linear()))
						rp_trial = rand() * r_LC
						p_max = maximum(dN_drp)
						y = rand() * p_max
						p_max = maximum(dN_drp)
						if y < itp_rp(rp_trial)
							push!(rp, rp_trial)
							push!(B, r_LC/rp_trial)
							break
						end
					end
				end
			end
			push!(accepted_rp,rp)
			push!(accepted_B,B)
		end
		push!(acepted_rp_total, accepted_rp)
		push!(acepted_B_total, accepted_B)
	end

	Specific_energies=vcat(vcat(accepted_energies_final...)...)*(params["sigma"]*(params["v_u"])/1000)^2 
	Beta=vcat(vcat(acepted_B_total...)...)

	a_km = - (GMsun_km * params["M"]) / (2 * Specific_energies)

	a_au = a_km / AU_in_km

	e=[]

	for i in 1:length(a_au)
		if a_au[i]<0
			e_i=1+(r_LC/Beta[i])/(abs(a_au[i]))
			push!(e,e_i)
		else
			e_i=1-(r_LC/Beta[i])/(a_au[i])
			push!(e,e_i)
		end
	end

	filename = "$(params["M"])_0.1AU_" * @sprintf("t1_%.3f_t2_%.3f_data_summary.txt", t1, t2)

	open(filename, "w") do io
		@printf(io, "Number of events between t = %.4f and %.4f Gyr: %d\n", t1, t2, N_events)
		@printf(io, "Method used: Explore snapshot by snapshot \n")
		@printf(io, "Mass of the MBH: %.2e solar masses\n\n", params["M"])
		
		@printf(io, "%-25s %-25s %-25s %-15s %-20s\n", "Specific energy [(km/s)^2]", "Penetration factor (Beta)", "Eccentricity", "Snapshot", "Snapshot time [Gyr]")
		@printf(io, "%s\n", "-"^115)

		for (e_val, b, ecc, snap_val, t_snap) in zip(Specific_energies, Beta, e, snap1, time_snap)
			@printf(io, "%-25.5e %-25.5e %-25.8f %-15d %-20.5f\n", e_val, b, ecc, snap_val, t_snap)
		end
	end
end

if method == 2
	function sample_interpolated_heatmap(Z, x, y, N)
		
		itp = interpolate((x, y), Z, Gridded(Linear()))
		
		function find_first_allzero_col_from_end(Z)
			ncols = size(Z, 2)
			for col_idx in ncols:-1:1
				if all(Z[:, col_idx] .== 0)
					return col_idx
				end
			end
			return 2  
		end

		idx = find_first_allzero_col_from_end(Z)

		xmin, xmax = minimum(x), maximum(x)
		ymin, ymax = y[idx-1], maximum(y)
		samples = []

		while length(samples) < N
			x₀ = rand() * (xmax - xmin) + xmin
			y₀ = rand() * (ymax - ymin) + ymin
			z = itp[x₀, y₀]
			if rand() < z / maximum(Z)
				push!(samples, (x₀, y₀))
			end
		end
		return samples
	end

	vals1=[]
	for i in 1:length(TDE_final[1])
		Evals2= sample_interpolated_heatmap(dNdrpdEdt[:,i,:], E, collect(r_p/r_LC), TDE_final[1][i])
		push!(vals1, Evals2 )
	end

	vals=vcat(vals1...)

	Evals = [-(x[1]-params["phis0"]) for x in vals]
	rpvals = [x[2] for x in vals]

	Specific_energies_2=Evals*(params["sigma"]*(params["v_u"])/1000)^2
	Beta_2=1 ./ rpvals

	a_km_2 = - (GMsun_km * params["M"]) / (2 * Specific_energies_2)

	a_au_2 = a_km_2 / AU_in_km

	e_2=[]

	for i in 1:length(a_au_2)
		if a_au_2[i]<0
			e_i=1+(r_LC/Beta_2[i])/(abs(a_au_2[i]))
			push!(e_2,e_i)
		else
			e_i=1-(r_LC/Beta_2[i])/(a_au_2[i])
			push!(e_2,e_i)
		end
	end


	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"Specific\ energy\ in\ internal\ units", ylabel=L"r_{p}/r_{lc}", limits=(-200,10,nothing,nothing), title="Distribution of E and rₚ between t = $t1 and t = $t2 Gyr")
	hm = heatmap!(ax, -(E.-params["phis0"]), r_p/r_LC, dNdrpdE; colormap=:viridis, colorrange=(0, 0.5*maximum(dNdrpdE)))
	Colorbar(fig[1,2], hm,  label = L"\frac{dN}{dE \, dr_p}")
	CairoMakie.scatter!(ax, Evals, rpvals, markersize=5, color=:red)

	save("E and rp dist with sampled points $(params["M"])_0.1AU.pdf", fig)

	filename = "$(params["M"])_0.1AU_" * @sprintf("t1_%.3f_t2_%.3f_data_summary.txt", t1, t2)

	open(filename, "w") do io
		@printf(io, "Number of events between t = %.4f and %.4f Gyr: %d\n", t1, t2, N_events)
		@printf(io, "Method used: Sampling from the heatmap integrate over the whole time-windows \n")
		@printf(io, "Mass of the MBH: %.2e solar masses\n\n", params["M"])
		
		@printf(io, "%-25s %-25s %-25s %-15s %-20s\n", "Specific energy [(km/s)^2]", "Penetration factor (Beta)", "Eccentricity", "Snapshot", "Snapshot time [Gyr]")
		@printf(io, "%s\n", "-"^115)

		for (e_val, b, ecc, snap_val, t_snap) in zip(Specific_energies_2, Beta_2, e_2, snap1, time_snap)
			@printf(io, "%-25.5e %-25.5e %-25.8f %-15d %-20.5f\n", e_val, b, ecc, snap_val, t_snap)
		end
	end
end