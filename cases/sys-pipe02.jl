include("../Dependency.jl")
using SYSApplication

input_file = "cases/sys-input/02single-pipe.json"
case = SYSApplication.Case("SYS-Pipe-02", input_file, "cases/output")

SYSApplication.initialize(case)
SYSApplication.apprun(case)
