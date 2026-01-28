First, the datasets are generated  (including both the dynamical systems and their corresponding labels) using .jl scripts for maps, flows, and colored noises.

Next, these generated time series are provided as input to the script Gera_Microestados_dados.jl, which computes the probabilities associated with each microstate considering microstates of size N = 2.

Finally, the systems themselves, and subsequently the microstate probability distributions, are used as inputs for the analyses performed in code.ipynb, where the MLP is built, trained and tested.
