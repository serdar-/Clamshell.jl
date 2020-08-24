module NetworkModels 

export GaussianNetworkModel, AnisotropicNetworkModel

abstract type AbstractNetworkModel end

struct GaussianNetworkModel <: AbstractNetworkModel
    r::Float64 # Radius (Å)
    Cα_coords::Array{Float64,2}
    eigvals::Array{Float64,1}
    eigvecs::Array{Float64,2}
    Laplacian::Array{Float64,2} 
end

struct AnisotropicNetworkModel <: AbstractNetworkModel
    γ::Float64 # Spring constant 
    r::Float64 # Radius (Å)
    Cα_coords::Array{Float64,2}
    eigvals::Array{Float64,1}
    eigvecs::Array{Float64,2}
    Hessian::Array{Float64,2} 
end

end