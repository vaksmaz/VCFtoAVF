from autopvs1 import AutoPVS1
demo=AutoPVS1('22-36678800-G-A')
print(demo.consequence, demo.hgvs_c, demo.hgvs_p, demo.pvs1.strength_raw,
      demo.pvs1.strength, demo.pvs1.criterion)
