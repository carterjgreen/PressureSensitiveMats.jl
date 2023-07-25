using PressureSensitiveMats
using Documenter

DocMeta.setdocmeta!(PressureSensitiveMats,
    :DocTestSetup,
    :(using PressureSensitiveMats);
    recursive = true)

makedocs(;
    modules = [PressureSensitiveMats],
    authors = "Carter Green <carter.green@icloud.com> and contributors",
    repo = "https://github.com/carterjgreen/PressureSensitiveMats.jl/blob/{commit}{path}#{line}",
    sitename = "PressureSensitiveMats.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://carterjgreen.github.io/PressureSensitiveMats.jl",
        assets = String[]),
    pages = ["Home" => "index.md",
        "combiners.md",
        "movement.md",
        "occupancy.md",
        "utils.md"])

deploydocs(;
    repo = "github.com/carterjgreen/PressureSensitiveMats.jl.git",
    devbranch = "main")
