using Documenter, formalA2

makedocs(sitename="My Documentation")
deploydocs(
    repo   = "https://github.com/mongibellili/formalA2.git",
    branch = "gh-pages",
)