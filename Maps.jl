using Random
using DelimitedFiles
using DynamicalSystems

Dir = "N2/Data"
mkpath(Dir)

rng          = MersenneTwister()
Random.seed!()

function main()

    Classes                               = 10    #Number of classes
    transient                             = 1000  
    data_size                             = 1000  #Time series size
    N_systems                             = 4000  
    
    Serie_In                              = rand(Float64,(N_systems,data_size)
    Labels                                = zeros(Int64,N_systems) 
    Vector_Labels                         = zeros(Int64,N_systems) 

    count                                 = 1
    for j=1:N_systems
        Labels[j]                         = count
        if (j > round(Int64,(count*N_systems/Classes)))
            count                         = count + 1
        end
    end

    for j=1:N_systems
        
        Int_Rand                          = rand(1:(N_systems-j+1))
        Label_Use                         = Labels[Int_Rand]
        Vector_Labels[j]                  = Label_Use

        deleteat!(Labels,Int_Rand)

        #Parameter Changes for BetaX
        β = 1.99 + (Label_Use - 1) * 0.5
        #Parameter Changes for LogisticMap
        #r                                 = 3.95+(Label_Use*(0.05/(1.0*Classes)))

        Out = 1.0
        println(j,' ',size(Labels),' ',r,' ',Label_Use)
        
        while (Out == 1.0)

            ####################DynamicalSystem(BetaX)##############################
            x                                 = rand()
            for i=1:(transient+data_size)
                x                             = β*x
                while (x > 1.00)
                    x                         = x-1.00
                end
                if (i > transient)
                    Serie_In[j,(i-transient)] = x
                end
            end

            if ((x > 0.0) && (x < 1.0))
                Out = 0.0
            end 
            ########################################################################
            ####################DynamicalSystem(Logistic Map)#######################
            #=x                                 = rand()
            for i=1:(transient+data_size)
                x                             = ((r*x)*(1.0-x))
                if (i > transient)
                    Serie_In[j,(i-transient)] = x
                end
            end

            if ((x > 0.0) && (x < 1.0))
                Out = 0.0
            end=#
            ########################################################################
            
        end
    end
    
    writedlm("$Dir/A_LabelsData_BetaX_10_Classes.dat",Vector_Labels)
    writedlm("$Dir/A_Data_BetaX_10_Classes.dat",Serie_In)

    #writedlm("$Dir/A_LabelsData_Logistic_10_Classes.dat",Vector_Labels)
    #writedlm("$Dir/A_Data_Logistic_10_Classes.dat",Serie_In)
    
end

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

return main()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
