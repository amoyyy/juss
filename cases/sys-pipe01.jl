include("../Dependency.jl")
using SYSApplication

input_file = "cases/sys-input/01single-pipe.json"
case = SYSApplication.Case("SYS-Pipe-01", input_file, "cases/output")

SYSApplication.initialize(case)
SYSApplication.apprun(case)
