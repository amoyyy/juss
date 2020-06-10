include("../Dependency.jl")
using SYSApplication

input_file = "cases/sys-input/IHX-02.json"
# 40 node, full flowrate

case = SYSApplication.Case("SYS-IHX-02", input_file, "cases/output")

SYSApplication.initialize(case)
SYSApplication.apprun(case)
