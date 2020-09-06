using JSExpr: @js

CANVAS_ID::Int64 = 0

function create_view(pdb_string::String,js::String)::HTML{String}
    template = """
    <script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js""></script>
    <div class="pdb_data" style="display:none;">$pdb_string</div>
    <div id="container-$CANVAS_ID" class="mol-container" style="width:700px;height:500px"></div>
    <script>
        let element = \$("#container-$CANVAS_ID");
        let config = { backgroundColor: '#fefcf5' };
        let v = \$3Dmol.createViewer( element, config );
        v.addModel( data, "pdb" );
        $js
    </script>
    """
    return HTML(template)
end

function set_style(style::Dict; selection::Dict=Dict())::String
    # v.setStyle({},{cartoon:{color:"black"}});
    selection_js = @js selection
    style_js = @js style
    style_string = """
    v.setStyle($selection,$style_js);"""
    return style_string
end






