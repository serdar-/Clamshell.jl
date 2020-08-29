
using Pkg
Pkg.activate("../")
push!(LOAD_PATH,"../src")
using Documenter, Clamshell

makedocs(sitename="Clamshell")