using Test
using BioStructures
using Distances
using NearestNeighbors


function test_download_pdb()
    downloadpdb("4AKE")
    true
end

function read_pdb()
    struc = read("./4AKE.pdb",PDB)
    cords = collectatoms(struc["A"],calphaselector) .|> 
                       coords |> (x) -> hcat(x...)
    return cords
end

function pw_distance()
    cords = read_pdb()
    r = pairwise(Euclidean(),cords,dims=1)
    print(size(r))
    true
end

function find_neighbors()
    cords = read_pdb()
    print(size(cords))
    r = 10. 
    balltree = BallTree(cords,Euclidean())
    idx = inrange(balltree,cords,r,false)
    print(idx)
    true
end

function test_remove_pdb()
    rm("4AKE.pdb")
    true
end

@test find_neighbors()
