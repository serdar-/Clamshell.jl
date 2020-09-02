using NearestNeighbors: BallTree, inrange, Euclidean
using LinearAlgebra: eigen
using StatsBase: corspearman

include("utils.jl")

"""
    AbstractNetworkModel

Abstract type used for defining network models. 
"""
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

function find_neighbors(Cα_coords::Array{Float64,2}; radius::Float64=10.)::Array{Array{Int64,1},1}
    balltree = BallTree(Cα_coords, Euclidean())
    neighbor_indices = inrange(balltree, Cα_coords, radius, false)
    return neighbor_indices
end

function build_Laplacian_matrix(Cα_coords::Array{Float64,2}; radius::Float64=7.3)::Array{Float64,2}
    N = size(Cα_coords)[2]
    laplacian = zeros(N, N)
    neighbors = find_neighbors(Cα_coords;radius=radius)
    @inbounds for i = 1:N
        laplacian[i,neighbors[i]] .= -1
        laplacian[i,i] = length(neighbors[i]) - 1
    end
    return laplacian
end

function GNM(Cα_coords::Array{Float64,2}; radius::Float64=7.3)::GaussianNetworkModel
    laplacian = build_Laplacian_matrix(Cα_coords; radius=radius)
    λ,Q = eigen(laplacian)
    gnm = GaussianNetworkModel(radius, Cα_coords, λ, Q, laplacian)
    return gnm
end

function GNM(ps::ProteinStructure; radius::Float64=7.3)::GaussianNetworkModel
    Cα_coords = collectatoms(ps,calphaselector) |> get_coords
    gnm = GNM(Cα_coords; radius=radius)
    return gnm
end

# To get hinges:
# corspearman(gnm.eigvecs[:,2]*gnm.eigvecs[:,2]'./gnm.eigvals[2])

function get_Hij(p1::Array{Float64,1}, p2::Array{Float64,1}, γ::Float64)::Array{Float64,2}
    if all(p1 == p2)
        return zeros(3, 3)
    else
        v = p1 - p2
        Hij = -γ .* (v * v') ./ (sum(v.^2))
        return Hij
    end
end

function build_Hessian_matrix(Cα_coords::Array{Float64,2}; radius::Float64=15.0,γ::Float64=1.)::Array{Float64,2}
    N = size(Cα_coords)[2]
    hessian = zeros(N * 3, N * 3)
    neighbors = find_neighbors(Cα_coords;radius=radius)
    ind = reshape([1:N * 3...], 3, N)
    @inbounds for i = 1:N 
        row = zeros(3, 3, N)
        @inbounds for j in neighbors[i]
            row[:,:,j] = get_Hij(Cα_coords[:,i], Cα_coords[:,j], γ)
        end
        hessian[ind[:,i],:] = reshape(row, 3, N * 3)
        hessian[ind[:,i],ind[:,i]] = -sum(row, dims=3) 
    end
    return hessian
end

function ANM(Cα_coords::Array{Float64,2}; radius::Float64=15.0, γ::Float64=1.)::AnisotropicNetworkModel
    hessian = build_Hessian_matrix(Cα_coords; radius=radius, γ=γ)
    λ,Q = eigen(hessian)
    anm = AnisotropicNetworkModel(γ,radius,Cα_coords,λ,Q,hessian)
    return anm
end

function ANM(ps::ProteinStructure; radius::Float64=15.0, γ::Float64=1.)::AnisotropicNetworkModel
    Cα_coords = collectatoms(ps,calphaselector) |> get_coords
    anm = ANM(Cα_coords; radius=radius, γ=γ)
    return anm
end

"""
    eigvals(nm::AbstractNetworkModel; n_modes::Int64=20)

Returns the eigenvalues calculated from the relevant elastic network model. 
"""
function eigvals(nm::AbstractNetworkModel; n_modes::Int64=20)::Array{Float64,1}
    if isa(nm,GaussianNetworkModel)
        if n_modes == -1
            return nm.eigvals[2:end]
        else
            return nm.eigvals[2:n_modes+1]
        end
    else
        if n_modes == -1
            return nm.eigvals[7:end]
        else
            return nm.eigvals[7:n_modes+6]
        end
    end
end

function eigvecs(nm::AbstractNetworkModel; n_modes::Int64=20)::Array{Float64,2}
    if isa(nm,GaussianNetworkModel)
        if n_modes == -1
            return nm.eigvecs[:,2:end]
        else
            return nm.eigvecs[:,2:n_modes+1]
        end
    else
        if n_modes == -1
            return nm.eigvecs[:,7:end]
        else
            return nm.eigvals[:,7:n_modes+6]
        end
    end
end

function mode_correlations(gnm::GaussianNetworkModel, mode::Int64=1)::Array{Float64,2}
    q = gnm.eigvecs[:,mode+1]
    λ = gnm.eigvals[mode+1]
    correlation = corspearman(q*q'./λ)
    return correlation
end