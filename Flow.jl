using Random
using DelimitedFiles
using DynamicalSystems
using DifferentialEquations
using Statistics
using FFTW
#using PyPlot
#####################################################################################################
############################################### Variables ############################################
#####################################################################################################

rng = MersenneTwister()
Random.seed!(rng)

function lorenz!(du, u, p, t)
    σ, ρ, β = p
    du[1] = σ * (u[2] - u[1])
    du[2] = u[1] * (ρ - u[3]) - u[2]
    du[3] = u[1] * u[2] - β * u[3]
end

function main(flow)
    # Parameters
    Classes       = 10
    transient     = 1000
    data_size     = 1000
    N_systems     = 4000
    runs          = 1 

    if flow == "Lorenz"
        system_func = lorenz!
        param_range = (σ = 10.0, ρ_base = 27.99, β = 8/3)  

    else
        error("Error")
    end

    for run in runs:runs
        Serie_In      = rand(Float64, (N_systems, data_size))  
        Labels        = zeros(Int64, N_systems)
        Vector_Labels = zeros(Int64, N_systems)

        # Labels
        count = 1
        for j in 1:N_systems
            Labels[j] = count
            if j > round(Int64, (count * N_systems / Classes))
                count += 1
            end
        end

        for j in 1:N_systems
            Int_Rand     = rand(1:(N_systems - j + 1))
            Label_Use    = Labels[Int_Rand]
            Vector_Labels[j] = Label_Use

            deleteat!(Labels, Int_Rand)

            if flow == "Lorenz"
                σ = param_range.σ
                ρ = param_range.ρ_base + Label_Use * (10.0 / (1.0*Classes))
                β = param_range.β
                parameters = [σ, ρ, β]
                u0 = [rand(), rand(), 0.0]
                h_step = 0.25

            end

            println("Run $run - Flow $flow: System $j, Label = $Label_Use")

            # Solve the system
            tspan = (0.0, ((data_size + transient) * h_step) - h_step)
            prob = ODEProblem(system_func, u0, tspan, parameters)
            sol = solve(prob, saveat=h_step)

	    # Save the x coordinate
            Serie_In[j,:]   = sol[1, (transient + 1):(transient + data_size)]
           end

        folder = "Data"
        isdir(folder) || mkpath(folder)

        # Salvar systems and labels
        writedlm("$folder/A_Data_$(flow)_10_Classes.dat", Serie_In)
        writedlm("$folder/Labels_$(flow)_10_Classes.dat", Vector_Labels)
    end
end


