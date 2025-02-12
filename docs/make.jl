using Documenter, DocumenterInterLinks

using Pkg
Pkg.add(url="https://github.com/ChristopheBerthod/Piecewise.jl", subdir="Piecewise")
Pkg.add(url="https://github.com/ChristopheBerthod/Piecewise.jl", subdir="PiecewiseHilbert")
Pkg.add(url="https://github.com/ChristopheBerthod/Piecewise.jl", subdir="PiecewiseLorentz")
using Piecewise, PiecewiseHilbert, PiecewiseLorentz

using MagnetoTransport

links = InterLinks(
    "Piecewise" => "https://ChristopheBerthod.github.io/Piecewise.jl/dev/objects.inv"
)

makedocs(
    repo = Documenter.Remotes.GitHub("ChristopheBerthod", "MagnetoTransport.jl"),
    sitename = "MagnetoTransport",
    format = Documenter.HTML(prettyurls = false,
        edit_link = "main",
        inventory_version = "v0.1.0"
        ),
    plugins = [links],
    modules = [MagnetoTransport],
    pages = [
        "index.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/ChristopheBerthod/MagnetoTransport.jl.git",
    devbranch = "main"
)
