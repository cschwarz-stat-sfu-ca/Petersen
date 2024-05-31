library(Matrix)
source("integrated_photoID_functions.r")

#----------------------------------------------------------------------------------------------------------------------#
#Example using the integrated photo-ID functions
#----------------------------------------------------------------------------------------------------------------------#

#Generate data 'Enc.Mat'
set.seed(34324)

N = 100           # population size
noccas = 2        # number of sampling occasions
p = 0.4           # probability of detection
delta_R = 0.4     # probability of right-sided encounter, given detection
delta_L = 0.4     # probability of left-sided encounter, given detection
alpha = 0.5       # probability of both-sided detection, given both sides were encountered (alpha=0 for LR data type and alpha=1 for LRB data type; this is enforced within the 'fun_sim_data' function no matter what is specified for alpha here)
alpha = 0
datatype = "LR"   # "LR" is LR data type; "LRB" is LRB data type; "LRAB" is LRAB data type

Enc.Mat=fun_sim_data(N,noccas,p,delta_R,delta_L,alpha,datatype)
Enc.Mat

# Get initial frequency count from the feasible set of latent history frequencies (x_0 in paper)
Freqs=fun_get_freq(Enc.Mat,datatype)
Freqs

# Get A matrix and latent histories
A_hists=fun_get_A(noccas,datatype)
A_hists

# Get basis vectors based on A and Freqs
Basis=fun_get_basis_vectors(noccas,A_hists$A,A_hists$ivect,Freqs,data.type=datatype)
Basis
