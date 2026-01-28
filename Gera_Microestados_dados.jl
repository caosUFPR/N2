using Random
using DelimitedFiles
using DynamicalSystems
using LinearAlgebra
using DataStructures
using Distributions
using Statistics
using SignalAnalysis
#using StatsBase

#####################################################################################################################################
Dir = "N2/Data"
mkpath(Dir)

rng          = MersenneTwister()
Random.seed!()

#####################################################################################################################################
###########################################################Variables S_Max###########################################################
#####################################################################################################################################

const StatsBlock         = 2 #Microstate size (this code works well until N=2, N=3 and N=4; for larger values, the algorithm still holds, but generally occurs problems with memory and enough data)
const Window_Size        = 1000 #Window size extracted from the data that will be evaluated
const Frac               = 10 #Split the interval of threshold (init_eps - max_eps) in Frac parts
const Frac2              = 5 #Given the interval where the maximum entropy is found, split in Frac2 parts to improve precision
const init_eps           = 0.00000001 #Minimum possible threshold
const max_eps            = 0.49999999 #Maximum possible threshold
const Max_Micro          = Int64(2^(StatsBlock*StatsBlock)) #Total number of microstates

const Sample_N           = floor(Int64,0.1*Window_Size*Window_Size) #Number of microstates samples
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

function Max_Entropy_Method(Serie,x_rand,y_rand,pow_vec)

    MicroStates             = zeros(Float64,Max_Micro)
    	
    S_Max=0.0; Threshold_Max=0.0 ; Threshold=init_eps; Var_Eps = (max_eps-init_eps)/Frac
    
    #Loop over different threholds (first part)
    for i=1:Frac2
        if (i > 1)
            Threshold=Threshold_Max-Var_Eps
            Var_Eps=2*Var_Eps/Frac
        end
	#Loop over different threholds (second part)
        for j=1:Frac
            Stats              	= zeros(Float64,Max_Micro);
	    #Loop that calculate the microstates probabilities for a given data and threshold
            for count=1:Sample_N
                Add=0
                for count_x=1:StatsBlock
                    for count_y=1:StatsBlock
                        if (abs(Serie[x_rand[count]+count_x]-Serie[y_rand[count]+count_y]) <= Threshold)
                            a_binary=1
                        else
                            a_binary=0
                        end
                        Add=Add+a_binary*pow_vec[count_y+((count_x-1)*StatsBlock)]
                    end
                end
                Stats[Int64(Add)+1]+=1
            end
            
            S=0
	    #Entropy Calculation
            for k=1:Max_Micro
                if (Stats[k] > 0)
                    S+=(-(Stats[k]/(1.0*Sample_N))*(log((Stats[k]/(1.0*Sample_N)))))
                end
            end
	    ####################
	    #Update of the maximum entropy and corresponding quantities
            if (S > S_Max)
                S_Max            = S
                Threshold_Max    = Threshold
                MicroStates[:]   = (Stats[:]./Sample_N)
            end
	    ####################
            Threshold=Threshold+Var_Eps #Scanning of different thresholds
        end
    end
    
    return S_Max,Threshold_Max,MicroStates
end

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

function main()
    
    x_rand            		= zeros(Int64,Sample_N)
    y_rand            		= zeros(Int64,Sample_N)
    
    Serie_In                  = readdlm("$Dir/A_Data_BetaX_10_Classes.dat")
                    
    N_systems                 = 4000
    times_size                = 1
    Jump                      = 1

    Classes                   = 10

                    
    MicroStates               = zeros(Float64,N_systems,times_size,Max_Micro); #Matrix with all microstates probabilities for each window
    S_Vec                     = zeros(Float64,N_systems,times_size) #Vector for maximum microstates entropy for each windowDa
    Eps_Vec                   = zeros(Float64,N_systems,times_size) #Vector for the threshold that reaches the maximum microstates entropy for each window

    Out_Vector                = zeros(Float64,N_systems,Max_Micro+2)

    for system_evaluation=1:N_systems
                 
        for times_size_change=1:times_size
                            
            Data_In                      = zeros(1,Window_Size)
            Data_In[1,:]                 = Serie_In[system_evaluation,(1+((times_size_change-1)*Window_Size)):Jump:((Jump*times_size_change*Window_Size))] #Window change in the time series
            Data_In[1,:]                 = ((Data_In[1,:].-minimum(Data_In))./(maximum(Data_In)-minimum(Data_In)))

            Data_Set_Analysis            = Dataset(transpose(Data_In))
                            
            x_rand = rand(1:(Window_Size - StatsBlock), Sample_N)
            y_rand = rand(1:(Window_Size - StatsBlock), Sample_N)
                    
            pow_vec = 2 .^ (0:(StatsBlock * StatsBlock - 1))
                            
            S_Vec[system_evaluation,times_size_change],Eps_Vec[system_evaluation,times_size_change],MicroStates[system_evaluation,times_size_change,:]   = Max_Entropy_Method(Data_In,x_rand,y_rand,pow_vec) #Calculation of the maximum recurrence microstates, threhold and microstates probabilities for each window

            Out_Vector[system_evaluation,1:Max_Micro]    = MicroStates[system_evaluation,times_size_change,:]
            Out_Vector[system_evaluation,Max_Micro+1]    = Eps_Vec[system_evaluation,times_size_change]
            Out_Vector[system_evaluation,Max_Micro+2]    = S_Vec[system_evaluation,times_size_change]
                            
            end

            println(system_evaluation)

        end
                    
        writedlm("$Dir/A_MicroStates_BetaX_10_Classes.dat",Out_Vector)            
    end

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

return main()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

