using Documenter, ModiaMath

makedocs(
  modules  = [ModiaMath],
  sitename = "ModiaMath",
  authors  = "Martin Otter (DLR-SR)",
  format = Documenter.HTML(prettyurls = false),
  pages    = [
     "Home"   => "index.md",
     "Manual" => [
        "man/Overview.md",
        "man/Plans.md"
        ],
     "Library" => [
        "lib/SimulationEngine.md"
        "lib/Result.md"
        "lib/DAE.md"
        "lib/NonlinearEquations.md"
        "lib/Variables.md"
        "lib/Logging.md"
        "lib/Frames.md"
        "lib/Utilities.md"
        "lib/ModiaMath.md"
     ]
  ]
)

