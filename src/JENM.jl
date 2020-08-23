module JENM

using NearestNeighbors
using BioStructures
using PROPACK
using LinearAlgebra

function get_calpha_atoms(ps::ProteinStructure;chain::String="A",model::Int64=1)::Array{AbstractAtom,1}
    cαs = collectatoms(ps[model][chain], calphaselector)
    return cαs
end

function get_calpha_atoms(ps::ProteinStructure; chains::Array{String,1}=["A"],models::Array{Int64,1}=[1])::Array{AbstractAtom,1}
    cαs = Array{AbstractAtom,1}()
    for model in models
        for chain in chains
            push!(cαs, collectatoms(ps[model][chain], calphaselector)...)  
        end
    end 
    return cαs
end

function get_calpha_coords(cαs::Array{AbstractAtom,1})::Array{Float64,2}
    cα_coords = cαs .|> coords |> (x) -> hcat(x...)
end

function find_neighbors(cα_coords::Array{Float64,2}; radius::Float64=10.)::Array{Array{Int64,1},1}
    balltree = BallTree(cα_coords, Euclidean())
    neighbor_indices = inrange(balltree, cα_coords, radius, false)
    return neighbor_indices
end

function build_Laplacian_matrix(cα_coords::Array{Float64,2}; radius::Float64=7.3)::Array{Float64,2}
    N = size(cα_coords)[2]
    laplacian = zeros(N, N)
    neighbors = find_neighbors(cα_coords;radius=radius)
    @inbounds for i = 1:N
        laplacian[i,neighbors[i]] .= -1
        laplacian[i,i] = length(neighbors[i]) - 1
    end
    return laplacian
end

function get_Hij(p1::Array{Float64,1},p2::Array{Float64,1},γ::Float64)::Array{Float64,2}
    if all(p1 == p2)
        return zeros(3,3)
    else
        v = p1 - p2
        Hij = -γ.*(v*v')./(sum(v.^2))
        return Hij
    end
end

function build_Hessian_matrix(cα_coords::Array{Float64,2}; radius::Float64=15.0,γ::Float64=1.)::Array{Float64,2}
    N = size(cα_coords)[2]
    hessian = zeros(N * 3, N * 3)
    neighbors = find_neighbors(cα_coords;radius=radius)
    ind = reshape([1:N*3...],3,N)
    @inbounds for i = 1:N 
        row = zeros(3,3,N)
        @inbounds for j in neighbors[i]
            row[:,:,j] = get_Hij(cα_coords[:,i],cα_coords[:,j],γ)
        end
        hessian[ind[:,i],:] = reshape(row,3,N*3)
        hessian[ind[:,i],ind[:,i]] = -sum(row,dims=3) 
    end
    return hessian
end

function decompose(c_matrix::Array{Float64,2})
    u, s, v = svd(c_matrix)
    return u, s, v
end

function decompose(c_matrix::Array{Float64,2}, n_modes::Int64)
    u, s, v, bnd, nprod, ntprod = tsvd_irl(c_matrix, k=n_modes)
    return u, s, v
end


end # module
