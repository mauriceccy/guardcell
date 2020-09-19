import scobra

phases = ['_Night', '_Open', '_Day', '_Close']
transitions = ['_NightOpen', '_OpenDay', '_DayClose', '_CloseNight']
frees = ['H2O_tx', 'CO2_tx', 'O2_tx', 'K_tx', 'Cl_tx'] # no H_tx

open_h = 1
day_h = 11
close_h = 1
night_h = 11

### osmolarity coefficients: KCl: ~0.9 (0.92); K2malate: ~0.8; sucrose ~1 (1.02)
oc_kcl = 0.9
oc_k2malate = 0.8
oc_sucrose = 1
osmol = 0.9 

m = scobra.Model("GC.xls")

for reac in m.Reactions("_biomass"):
	m.SetConstraint(reac, 0, 0)
for reac in m.Reactions("_tx"):
	m.SetConstraint(reac, 0, 0)

photon = [0, None]
m.SetConstraint("Photon_tx_Open", photon[0], photon[1])
m.SetConstraint("Photon_tx_Day", photon[0], photon[1])
m.SetReacsFixedRatio({"Photon_tx_Open":open_h, "Photon_tx_Day":day_h})

for p in phases:
	for f in frees:
		m.SetConstraint(f+p, None, None)
	m.SetConstraint("NADPH_Dehydrogenase_p"+p, 0, 0)
	m.SetConstraint("Plastoquinol_Oxidase_p"+p, 0, 0)
	m.SetReacsFixedRatio({"RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"+p:3, "RXN_961_p"+p:1})
	#m.SetObjective('Photon_tx' + p)

for t in transitions:
	m.SetSumReacsConstraint({"Cl"+t+"_storage":-1, "Malate"+t+"_storage":-2, "K"+t+"_storage":1}, 0)

m.SetSumReacsConstraint({"Cl_NightOpen_storage":2*oc_kcl, "Malate_NightOpen_storage":3*oc_k2malate, "Sucrose_NightOpen_storage":1*oc_sucrose}, 0)
m.SetSumReacsConstraint({"Cl_OpenDay_storage":2*oc_kcl, "Malate_OpenDay_storage":3*oc_k2malate, "Sucrose_OpenDay_storage":1*oc_sucrose}, osmol)
m.SetSumReacsConstraint({"Cl_DayClose_storage":2*oc_kcl, "Malate_DayClose_storage":3*oc_k2malate, "Sucrose_DayClose_storage":1*oc_sucrose}, osmol)
m.SetSumReacsConstraint({"Cl_CloseNight_storage":2*oc_kcl, "Malate_CloseNight_storage":3*oc_k2malate, "Sucrose_CloseNight_storage":1*oc_sucrose}, 0)

#restrict cl
m.SetConstraint('Cl_OpenDay_storage', 0, 0)
m.SetConstraint('Cl_DayClose_storage', 0, 0)
m.SetConstraint('Cl_CloseNight_storage', 0, 0)
m.SetConstraint('Cl_NightOpen_storage', 0, 0)

#m.SetConstraint("K_NightOpen_storage", 0 ,0)
#m.SetConstraint("K_OpenDay_storage", 0, 0)
m.SetConstraint("K_DayClose_storage", 0, 0)
#m.SetConstraint("K_CloseNight_storage", 0, 0)

#m.SetConstraint("Sucrose_NightOpen_storage", 0, 0)
m.SetConstraint("Sucrose_OpenDay_storage", 0, 0)
#m.SetConstraint("Sucrose_DayClose_storage", 0, 0)
#m.SetConstraint("Sucrose_CloseNight_storage", 0, 0)

#restrict starch
#m.SetConstraint('Starch_OpenDay_storage', 0, 0)
#m.SetConstraint('Starch_DayClose_storage', 0, 0)
#m.SetConstraint('Starch_CloseNight_storage', 0, 0)
#m.SetConstraint('Starch_NightOpen_storage', 0, 0)

#restrict malate
#m.SetConstraint('Malate_OpenDay_storage', 0, 0)
#m.SetConstraint('Malate_DayClose_storage', 0, 0)
#m.SetConstraint('Malate_CloseNight_storage', 0, 0)
#m.SetConstraint('Malate_NightOpen_storage', 0, 0)

atp = 2.866827556863427 # match carbon fixation = 20% of osmolyte accum
m.SetConstraint("ATPase_tx_Open", atp*open_h, atp*open_h)
m.SetConstraint("ATPase_tx_Day", atp*day_h, atp*day_h)
m.SetConstraint("ATPase_tx_Close", atp*close_h, atp*close_h)
m.SetConstraint("ATPase_tx_Night", atp*night_h, atp*night_h)
nadph = 9.0
for phase in phases:
	if phase == '_Night':
		m.SetConstraint('NADPHoxc_tx' + phase, (atp*night_h/nadph), (atp*night_h/nadph))
		m.SetConstraint('NADPHoxp_tx' + phase, (atp*night_h/nadph), (atp*night_h/nadph))
		m.SetConstraint('NADPHoxm_tx' + phase, (atp*night_h/nadph), (atp*night_h/nadph))
	if phase == '_Open':
		m.SetConstraint('NADPHoxc_tx' + phase, (atp*open_h/nadph), (atp*open_h/nadph))
		m.SetConstraint('NADPHoxp_tx' + phase, (atp*open_h/nadph), (atp*open_h/nadph))
		m.SetConstraint('NADPHoxm_tx' + phase, (atp*open_h/nadph), (atp*open_h/nadph))
	if phase == '_Day':
		m.SetConstraint('NADPHoxc_tx' + phase, (atp*day_h/nadph), (atp*day_h/nadph))
		m.SetConstraint('NADPHoxp_tx' + phase, (atp*day_h/nadph), (atp*day_h/nadph))
		m.SetConstraint('NADPHoxm_tx' + phase, (atp*day_h/nadph), (atp*day_h/nadph))
	if phase == '_Close':
		m.SetConstraint('NADPHoxc_tx' + phase, (atp*close_h/nadph), (atp*close_h/nadph))
		m.SetConstraint('NADPHoxp_tx' + phase, (atp*close_h/nadph), (atp*close_h/nadph))
		m.SetConstraint('NADPHoxm_tx' + phase, (atp*close_h/nadph), (atp*close_h/nadph))

sucrose = 0.72 # 80% of osmolyte from sucrose
m.SetConstraint('Sucrose_tx_Day', sucrose, sucrose)
m.MinFluxSolve()
sol = m.GetSol(IncZeroes=True)
m.PrintSol('storage')
m.PrintSol('Photon')
m.PrintSol('CO2_tx')
sol.AsMtx().ToFile('GC_standard_minflux.csv')

m.SetBounds(1000) # to set a finite bound to bypass cobrapy unbounded problem
fva = m.FVA()
df = fva.AsMtx().sort_index()
df.to_csv('GC_standard_minflux_FVA.csv')
