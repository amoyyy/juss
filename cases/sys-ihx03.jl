include("../Dependency.jl")
using SYSApplication

input_file = "cases/sys-input/IHX-03.json"
# 40 node, full flowrate

case = SYSApplication.Case("SYS-IHX-03", input_file, "cases/output")

SYSApplication.initialize(case)
SYSApplication.apprun(case)
