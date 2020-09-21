using Base: IOBuffer
using Printf
using JSExpr: @js_str
using BioStructures: ProteinStructure, 
                     collectatoms, 
                     writepdb,
                     AbstractAtom,
                     Atom
using StatsBase: mean

include("network_models.jl")

CANVAS_ID = 0

function create_structure_view(pdb_string::String,js::String)::HTML{String}
    template = """
    <script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js""></script>
    <div class="pdb_data_$CANVAS_ID" style="display:none;">$pdb_string</div>
    <div id="container-$CANVAS_ID" class="mol-container" style="width:700px;height:500px"></div>
    <script>
        let element = \$("#container-$CANVAS_ID");
        let config = { backgroundColor: '#fefcf5' };
        let v = \$3Dmol.createViewer( element, config );
        let data = \$(".pdb_data_$CANVAS_ID").html();
        let model;
        if(data.length != 0){
            model = v.addModel( data, "pdb" );
        }
        $js
        v.zoomTo(); 
        v.render();
        v.zoom(1.1, 750);
    </script>
    """
    return HTML(template)
end

function create_network_view(js::String)::HTML{String}
    template = """
    <script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js""></script>
    <div id="container-$CANVAS_ID" class="mol-container" style="width:700px;height:500px"></div>
    <script>
        let element = \$("#container-$CANVAS_ID");
        let config = { backgroundColor: '#fefcf5' };
        let v = \$3Dmol.createViewer( element, config );
        $js
    </script>
    """
    return HTML(template)
end

function set_style(style::Dict; selection::Dict=Dict())::String
    # v.setStyle({},{cartoon:{color:"black"}});
    selection_js = js"$selection"
    style_js = js"$style"
    style_string = """
    model.setStyle($selection_js,$style_js);"""
    return style_string
end

function create_pdb_string(atoms::Array{AbstractAtom,1})::String
    io = IOBuffer()
    writepdb(io,atoms)
    pdb_string = String(take!(io))
    return pdb_string
end

function show_structure(ps::ProteinStructure; 
                        model::Int64=1,
                        chains::Array{String}=["A"],
                        style::Dict=Dict("cartoon"=>Dict("color"=>"#5e7ad3")))::HTML{String}
    pdb_string =  ps[model][chains...] |> collectatoms |> create_pdb_string
    style_string = set_style(style)
    global CANVAS_ID += 1
    view = create_structure_view(pdb_string,style_string)
    return view
end
# Cα_coords
function show_correlations(atoms::Array{AbstractAtom,1};mode::Int64=1,show_hinges::Bool=false,additional_style::String="")::HTML{String}
    Cα_coords = get_calpha_coords(atoms)
    gnm = GNM(Cα_coords) 
    corrs = mode_correlations(gnm,mode) |> (x) -> x[:,1]
    cₚ = findall(round.(corrs) .== 1.0) 
    cₙ = findall(round.(corrs) .== -1.0)
    cₚ_style = Dict("cartoon"=>Dict("color"=>"#5e7ad3"))
    cₙ_style = Dict("cartoon"=>Dict("color"=>"#f57464"))
    style_string = set_style(cₚ_style; selection=Dict("resi"=>cₚ))
    style_string *= set_style(cₙ_style; selection=Dict("resi"=>cₙ))
    if show_hinges
        hinges = get_hinge_indices(gnm)
        hinge_style = Dict("cartoon"=>Dict("color"=>"#ce7d0a"))
        style_string *= set_style(hinge_style;
                                  selection=Dict("resi"=>hinges))
    end
    pdb_string = create_pdb_string(atoms)
    global CANVAS_ID += 1
    structure_view = create_structure_view(pdb_string,style_string*additional_style)
    return structure_view
end

function convert_coords(x::Array{Float64,2},i::Int64)::Dict
    return Dict("x"=>x[1,i],"y"=>x[2,i],"z"=>x[3,i])
end

function convert_coords(x::Array{Float64,2})::Dict
    return Dict("x"=>x[1,1],"y"=>x[2,1],"z"=>x[3,1])
end

function show_network(atoms::Array{AbstractAtom,1};radius::Float64=7.3,show_structure::Bool=false)::HTML{String}
    Cα_coords = get_calpha_coords(atoms)
    neighbors = find_neighbors(Cα_coords;radius=radius);
    if show_structure
        pdb_string = create_pdb_string(atoms)
        model_style = set_style(Dict("cartoon"=>Dict("color"=>"#5e7ad3")))
    else
        pdb_string = ""
        model_style = ""
    end
    shape = "let network = v.addShape(1, {color:'red'});"
    N = size(Cα_coords)[2]
    # convert_coords = (x,i) -> Dict("x"=>x[1,i],"y"=>x[2,i],"z"=>x[3,i])
    @inbounds for i = 1:N
        @inbounds for j = neighbors[i]
            if i != j
                p_start = convert_coords(Cα_coords,i) |> (x) -> @sprintf("%s",js"$x")
                p_end = convert_coords(Cα_coords,j) |> (x) -> @sprintf("%s",js"$x")
                shape *= "network.addCylinder({start:$p_start,end:$p_end,radius:0.1,fromCap:2,toCap:2,color:'red',dashed:false});"
            end
        end
    end
    global CANVAS_ID += 1
    network_view = create_structure_view(pdb_string,shape*model_style)
    return network_view
end

function show_hinge_plane(atoms::Array{AbstractAtom,1};radius::Float64=7.3,show_structure::Bool=true)::HTML{String}
    Cα_coords = get_calpha_coords(atoms)
    gnm = GNM(Cα_coords;radius=radius)
    n_vector = hinge_plane_normal(gnm) # Normal of hinge plane
    hinges = get_hinge_indices(gnm)
    hinges_mean = mean(Cα_coords[:,hinges],dims=2)
    if show_structure
        pdb_string = create_pdb_string(atoms)
        model_style = set_style(Dict("cartoon"=>Dict("color"=>"#5e7ad3")))
    else
        pdb_string = ""
        model_style = ""
    end
    p_start = hinges_mean |> convert_coords |> (x) -> @sprintf("%s",js"$x")
    p_end = hinges_mean + n_vector.*10 |> convert_coords |> (x) -> @sprintf("%s",js"$x")
    plane = """let shape = v.addShape({color:'#8c65bd',alpha:0.7});
    let plane = shape.addCylinder({start:$p_start,end:$p_end,fromCap:1,toCap:1,radius:20.0});"""
    global CANVAS_ID += 1
    # show_correlations(atoms::Array{AbstractAtom,1};mode::Int64=1,show_hinges::Bool=false,additional_style::String="")
    correlation_view = show_correlations(atoms;additional_style=plane)
    return correlation_view
end