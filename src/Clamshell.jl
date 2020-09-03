module Clamshell

export AbstractNetworkModel,
       AnisotropicNetworkModel,
       GaussianNetworkModel,
       GNM,
       ANM,
       eigvals,
       eigvecs,
       get_calpha_coords,
       mode_correlations

include("network_models.jl")

end