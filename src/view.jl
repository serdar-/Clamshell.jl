using JSExpr
using BioStructures

function read_pdb()
    ps = read("./4AKE.pdb",PDB)
    return ps 
end

function convert_atom_to_js(atom)
    atom_dict = Dict(
                    :resn => strip(atom.residue.name),
                    :resi => atom.residue.number,
                    :x => atom.coords[1],
                    :y => atom.coords[2],
                    :z => atom.coords[3],
                    :elem => strip(atom.element)
                    )
    return @js $atom_dict
end

function test_view()
    HTML(@js """<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.5.4/3Dmol.js" integrity="sha512-u0f/7wxqD3UY79xCVst5bxqPDzzXJcXDvsIaM39hvhhuXaTAMBwT2izHm3ynSNSgH3MvTkT5CEIaTWQaa/Sh6g==" crossorigin="anonymous"></script>
    <div id="container-01" class="mol-container" style="width:540px;height:540px;"></div>
    <script>
        let viewer;
        let element = \$('#container-01');
        let config = { backgroundColor: 'white' };
        viewer = \$3Dmol.createViewer( element, config );
        viewer.addSphere({ center: {x:0, y:0, z:0}, radius: 10.0, color: 'green' });
        viewer.zoomTo();
        viewer.render();
        viewer.zoom(0.8, 2000);
    </script>""")
end
#     /**
#  * Atom representation. Depending on the input file format, not all fields may be defined.
#  * @typedef AtomSpec
#  * @prop {string} resn - Parent residue name
#  * @prop {number} x - Atom's x coordinate
#  * @prop {number} y - Atom's y coordinate
#  * @prop {number} z - Atom's z coordinate
#  * @prop {ColorSpec} color - Atom's color, as hex code or built-in color string
#  * @prop {ColorSpec} surfaceColor - Hex code for color to be used for surface patch over this atom
#  * @prop {string} elem - Element abbreviation (e.g. 'H', 'Ca', etc)
#  * @prop {boolean} hetflag - Set to true if atom is a heteroatom
#  * @prop {string} chain - Chain this atom belongs to, if specified in input file (e.g 'A' for chain A)
#  * @prop {number} resi - Residue number
#  * @prop {number} icode
#  * @prop {number} rescode
#  * @prop {number} serial - Atom's serial id number
#  * @prop {string} atom - Atom name; may be more specific than 'elem' (e.g 'CA' for alpha carbon)
#  * @prop {Array.<number>} bonds - Array of atom ids this atom is bonded to
#  * @prop {string} ss - Secondary structure identifier (for cartoon render; e.g. 'h' for helix)
#  * @prop {boolean} singleBonds - true if this atom forms only single bonds or no bonds at all
#  * @prop {Array.<number>} bondOrder - Array of this atom's bond orders, corresponding to bonds identfied by 'bonds'
#  * @prop {Object} properties - Optional mapping of additional properties
#  * @prop {number} b - Atom b factor data
#  * @prop {string} pdbline - If applicable, this atom's record entry from the input PDB file (used to output new PDB from models)
#  * @prop {boolean} clickable - Set this flag to true to enable click selection handling for this atom
#  * @prop {function(this, $3Dmol.GLViewer)} callback - Callback click handler function to be executed on this atom and its parent viewer
#  * @prop {boolean} invert - for selection, inverts the meaning of the selection
#  */
end