# julia --color=yes --project make.jl
# julia -e 'using LiveServer; serve(dir="build")'
#
# refs: 
# [1] https://documenter.juliadocs.org/stable/man/guide/#Package-Guide
# [2] https://docs.julialang.org/en/v1/manual/documentation/
# [3] https://juliadocs.org/DocumenterCitations.jl/stable/

push!(LOAD_PATH, "../src/");

using Documenter, DocumenterCitations, MembraneAleFem, Pkg


PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION      = PROJECT_TOML["version"]
NAME         = PROJECT_TOML["name"]
AUTHORS      = join(PROJECT_TOML["authors"], ", ")
GITHUB       = "https://github.com/sahu-lab/MembraneAleFem.jl"

Manual = [
  "man/overview.md",
  "man/input.md",
  "man/analysis.md",
  "man/output.md"
]

const PAGES = [
  "Home" => "index.md",
  "Manual" => Manual,
]

bib = CitationBibliography(
  joinpath(@__DIR__, "src", "refs.bib"),
  style=:authoryear,
)

makedocs(
  sitename  = "$NAME.jl",
  authors   = "$AUTHORS",
  pages     = PAGES,
  format    = Documenter.HTML(
    prettyurls = true,
    canonical  = "https://sahu-lab.github.io/MembraneAleFem.jl",
    assets     = ["assets/favicon.ico"],
    footer     = "[$NAME.jl]($GITHUB) v$VERSION docs powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl)."
  ),
  plugins   = [bib],
);

deploydocs(
  repo="github.com/sahu-lab/MembraneAleFem.jl.git",
);

