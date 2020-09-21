using Pkg
Pkg.activate("./")
push!(LOAD_PATH,"../src")
using Documenter, Clamshell

makedocs(
    format=Documenter.HTML(assets=[asset("https://3Dmol.org/build/3Dmol-min.js")],prettyurls=false),
    sitename="Clamshell",
    pages = [
        "Home" => "index.md",
        "Introduction" => [
            "Brief introduction to ENM" => "intro_to_enm.md"
        ]
    ]
)
