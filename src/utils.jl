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

function get_coords(atoms::Array{AbstractAtom,1})::Array{Float64,2}
    atom_coords = atoms .|> coords |> (x) -> hcat(x...)
end

function decompose(c_matrix::Array{Float64,2})
    u, s, v = svd(c_matrix)
    return u, s, v
end

function decompose(c_matrix::Array{Float64,2}, n_modes::Int64)
    u, s, v, bnd, nprod, ntprod = tsvd_irl(c_matrix, k=n_modes)
    return u, s, v
end

