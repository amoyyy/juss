include("../Dependency.jl")
using SYSApplication

input_file = "cases/sys-input/IHX-01.json"
# 20 node, full flowrate

case = SYSApplication.Case("SYS-IHX-01", input_file, "cases/output")

#SYSApplication.initialize(case)
#SYSApplication.apprun(case)
code_native(SYSApplication.initialize)
