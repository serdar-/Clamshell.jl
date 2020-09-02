using BioStructures: ProteinStructure, calphaselector, collectatoms, coords

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

function get_calpha_coords(ps::ProteinStructure; chains::Array{String,1}=["A"],models::Array{Int64,1}=[1])::Array{Float64,2}
    return get_calpha_atoms(ps;chains=chains,models=models) |> get_coords
end
