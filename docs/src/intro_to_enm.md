In elastic network models, protein structures are simplified by considering only the α-carbon atoms of the protein structure. Thus each residue is represented by its α-carbon atom. Let's look at the example of pyruvate kinase structure. α-carbon atoms are shown in black spheres below.


```@raw html
<div id="container-01" class="mol-container" style="width:540px;height:540px;position:relative;"></div>
<script>
    let viewer;
    let element = $('#container-01');
    let config = { backgroundColor: 'white' };
    viewer = $3Dmol.createViewer( element, config );
    $3Dmol.download("pdb:1PKL",
                    viewer,{},
                    function(){
                        let model = viewer.getModel();
                        let atoms = viewer.selectedAtoms({chain:'A',invert:true});
                        model.removeAtoms(atoms);
                        //viewer.setStyle({chain:'A'},{cartoon:{color:'#5e7ad3',opacity:0.5}});
                        viewer.setStyle({chain:'A'},{cartoon:{color:'#5e7ad3',opacity:0.8}});
                        viewer.addStyle({chain:'A',atom:'CA'},{sphere:{radius:0.5,color:'black'}});
                        viewer.center({chain:'A',atom:'CA'});
                        viewer.zoomTo(1.1);
                        viewer.render();
                    });
</script>
```
These α-carbon atoms are going to interact with each other within a cutoff radius ($r_c$) that we define. As an example, let us consider that $r_c = 7.3Å$.

```@raw html
<div id="container-02" class="mol-container" style="width:540px;height:540px;position:relative;"></div>
<script>
    let v2;
    let e2 = $('#container-02');
    v2 = $3Dmol.createViewer( e2, config );
    distance = (a1,a2) => {
        dx = Math.pow(a1.x - a2.x,2);
        dy = Math.pow(a1.y - a2.y,2);
        dz = Math.pow(a1.z - a2.z,2);
        return Math.pow(dx + dy + dz,0.5);
    }
    $3Dmol.download("pdb:1PKL",
                    v2,{},
                    function(){
                        let model = v2.getModel();
                        let atoms = v2.selectedAtoms({chain:'A',invert:true});
                        let caAtoms = v2.selectedAtoms({chain:'A',atom:'CA'});
                        model.removeAtoms(atoms);
                        //viewer.setStyle({chain:'A'},{cartoon:{color:'#5e7ad3',opacity:0.5}});
                        for(let i; i < atoms.length; i++){
                            for(let j; j < i; j++){
                                var d = distance(atoms[i],atoms[j]);
                                if(d <= 7.3){
                                    v2.addCylinder({start:{x:atoms[i].x,y:atoms[i].y,z:atoms[i].z},
                                  end:{x:atoms[j].x,y:atoms[j].y,z:atoms[j].z},
                                  radius:0.5,
                                  fromCap:1,
                                  toCap:1});
                                }
                            }
                        }
                        v2.setStyle({chain:'A'},{cartoon:{color:'#5e7ad3',opacity:0.8}});
                        v2.addStyle({chain:'A',atom:'CA'},{sphere:{radius:0.5,color:'black'}});
                        v2.center({chain:'A',atom:'CA'});
                        v2.zoomTo(1.1);
                        v2.render();
                    });
</script>
```