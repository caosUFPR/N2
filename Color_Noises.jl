using Random
using DelimitedFiles
using FFTW 
include("NoiseGenerator.jl")

Dir = "N2/Data"
mkpath(Dir)

rng = MersenneTwister()
Random.seed!(rng)

function normalize(series)
    min_val = minimum(series)
    max_val = maximum(series)
    return (series .- min_val) ./ (max_val - min_val)
end

function main()
    Classes = 10
    transient = 10  
    data_size = 1000 
    N_systems = 4000  

    Serie_In = rand(Float64, (N_systems, data_size))
    Labels = zeros(Int64, N_systems)
    Vector_Labels = zeros(Int64, N_systems) 
    alphas = -2 .+ (0:Classes-1) * (4 / (Classes - 1))
    
    for j in 1:N_systems

        class_index = ((j - 1) % Classes) + 1
        alpha = alphas[class_index]
        Labels[j] = class_index
        Vector_Labels[j] = class_index

        x = NoiseGenerator(data_size, alpha)
        Serie_In[j, :] = normalize(x) 

        println(j, " ", length(Labels) - j, " ", alpha, " ", Labels[j])
    end


    perm = randperm(rng, N_systems)
    Serie_In = Serie_In[perm, :]
    Vector_Labels = Vector_Labels[perm]

    writedlm("$Dir/A_Data_Color_Noises_10_Classes.dat", Serie_In)
    writedlm("$Dir/A_LabelsData_Color_Noises_10_Classes.dat", Vector_Labels)

end

