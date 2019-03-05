from module_description_parser import ModuleDescription
s = "K00330+K00343"
module = ModuleDescription(s)
print(module.num_steps())
print(module.num_covered_steps(["K00330", "K00343"]))
print(module.num_covered_steps(["K00330"]))
print(module.num_covered_steps(["K00343"]))
module.kos()
print("Done")
