import numpy as np
import matplotlib.pyplot as plt
import flopy
import os
import xmf6

phys = dict(
    specific_discharge = 0.1,  # Specific discharge ($cm s^{-1}$)
    hydraulic_conductivity = 0.01,  # Hydraulic conductivity ($cm s^{-1}$)
    source_concentration = 1.0,  # Source concentration (unitless)
)
xmf6.nice_print(phys, "Parámetros físicos")

# Parámetros de la simulación (flopy.mf6.MFSimulation)
init = {
    'sim_name' : "flow",
#    'exe_name' : "C:\\Users\\luiggi\\Documents\\GitSites\\mf6_tutorial\\mf6\\windows\\mf6",
    'exe_name' : "../../mf6/macosarm/mf6",
    'sim_ws' : "output_flow_1D"
}

# Parámetros para el tiempo (flopy.mf6.ModflowTdis)
tdis = {
    'units': "seconds",
    'nper' : 1,
    'perioddata': [(120.0, 1, 1.0)]
}

# Parámetros para la solución numérica (flopy.mf6.ModflowIms)
ims = {}

# Parámetros para el modelo de flujo (flopy.mf6.ModflowGwf)
gwf = { 
    'modelname': init["sim_name"],
#    'model_nam_file': f"{init["sim_name"]}.nam",
    'save_flows': True
}

# Parámetros para la discretización espacial (flopy.mf6.ModflowGwfdis)
dis = {
    'length_units' : "centimeters",
    'nlay': 1, 
    'nrow': 1, 
    'ncol': 120,
    'delr': 0.1, 
    'delc': 0.1, 
    'top' : 1.0, 
    'botm': 0.0 
}

# Parámetros para las condiciones iniciales (flopy.mf6.ModflowGwfic)
ic = {
    'strt': 1.0
}

# Parámetros para las condiciones de frontera (flopy.mf6.ModflowGwfchd)
chd = {
    'stress_period_data': [[(0, 0, dis['ncol'] - 1), 1.0]],     
}

# Parámetros para las propiedades de flujo (flopy.mf6.ModflowGwfnpf)
npf = {
    'save_specific_discharge': True,
    'save_saturation': True,
    'icelltype':  0,
    'k': 0.01  # Hydraulic conductivity ($cm s^{-1}$) 
}

phys = dict(
    specific_discharge = 0.1,  # Specific discharge ($cm s^{-1}$)
    hydraulic_conductivity = 0.01,  # Hydraulic conductivity ($cm s^{-1}$)
    source_concentration = 1.0,  # Source concentration (unitless)
)
# Parámetros para las propiedades de los pozos (flopy.mf6.ModflowGwfwel)
phys["specific_discharge"] = 0.1 # Specific discharge ($cm s^{-1}$)
phys["source_concentration"] = 1.0  # Source concentration (unitless)
q = phys["specific_discharge"] * dis['delc'] * dis['delr'] * dis['top']

well = {
    'stress_period_data': [[(0, 0, 0), q, phys["source_concentration"],]],
    'pname': "WEL-1",
    'auxiliary' : ["CONCENTRATION"],
#    'save_flows': True
}

# Parámetros para almacenar y mostrar la salida de la simulación (flopy.mf6.ModflowGwfoc)
oc = {
    'budget_filerecord': f"{init['sim_name']}.bud",
    'head_filerecord': f"{init['sim_name']}.hds",
    'saverecord': [("HEAD", "ALL"), ("BUDGET", "ALL")],
#    'printrecord': [("HEAD", "ALL")]
}

# --- Inicialización de la simulación ---
o_sim = xmf6.common.init_sim(silent = True, init = init, tdis = tdis, ims = ims)
o_gwf, packages = xmf6.gwf.set_packages(o_sim, silent = True,
                                        gwf = gwf, 
                                        dis = dis, ic = ic, chd = chd, npf = npf, oc = oc, well = well)

# --- Escritura de archivos ---
o_sim.write_simulation(silent = True)

# --- Ejecución de la simulación ---
o_sim.run_simulation(silent = True)

# --- Recuperamos los resultados de la simulación ---
head = xmf6.gwf.get_head(o_gwf)
qx, qy, qz, n_q = xmf6.gwf.get_specific_discharge(o_gwf, text="DATA-SPDIS")

x, y, z = o_gwf.modelgrid.xyzcellcenters
row_length = o_gwf.modelgrid.extent[1]

# --- Definición de la figura ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize =(8,4), height_ratios=[2,0.1])

# --- Gráfica 1. ---
#plt.rcParams['font.family'] = 'DeJavu Sans'
ax1.plot(x[0], head[0, 0], marker=".", ls ="-", mec="blue", mfc="none", markersize="1", label = 'Head')
ax1.set_xlim(0, 12)
ax1.set_xticks(ticks=np.linspace(0, row_length,13))
ax1.set_xlabel("Distance (cm)")
ax1.set_ylabel("Head")
ax1.grid()

# --- Gráfica 2. ---
pmv = flopy.plot.PlotMapView(o_gwf, ax=ax2)
pmv.plot_grid(colors='dimgray', lw=0.5)
ax2.set_yticks(ticks=[0, 0.1])#, fontsize=8)
ax2.set_title("Mesh")
plt.tight_layout()
plt.show()