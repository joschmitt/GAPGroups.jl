using Documenter, GAPGroups

DocMeta.setdocmeta!( GAPGroups, :DocTestSetup, :( using GAPGroups ); recursive = true )

makedocs(
         format   = Documenter.HTML(),
         sitename = "GAPGroups.jl",
         modules = [GAPGroups],
         clean = true,
         doctest = true,
         strict = false,
         checkdocs = :none,
         pages    = [
             "index.md",
         ]
)

deploydocs(
   repo   = "github.com/oscar-system/GAPGroups.jl.git",
   target = "build",
)
