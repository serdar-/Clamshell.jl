module Clamshell

include("network_models.jl")

export AbstractNetworkModel,
       AnisotropicNetworkModel,
       GaussianNetworkModel,
       GNM,
       ANM,
       eigvals,
       eigvecs
end