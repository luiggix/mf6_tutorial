import os, sys
import numpy as np
import matplotlib.pyplot as plt
import flopy
from modflowapi import ModflowApi # La API
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
    'exe_name' : "C:\\Users\\luiggi\\Documents\\GitSites\\mf6_tutorial\\mf6\\windows\\mf6",
#    'exe_name' : "../../mf6/macosarm/mf6",
    'sim_ws' : "sandbox4"
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
#o_sim.run_simulation(silent = True)

###########

OS = sys.platform
if OS == "win32":
    mf6_lib = "libmf6.dll"
elif OS == "darwin":
    mf6_lib = "libmf6.dylib"
else:
    mf6_lib = "libmf6.so"

# Rutas a la biblioteca compartida y al archivo de configuración
mf6_lib_path = os.path.abspath(os.path.join("..", "..", "mf6", "windows", mf6_lib))
mf6_config_file = os.path.join(o_sim.sim_path, 'mfsim.nam')
print("Shared library:", mf6_lib_path)
print("Config file:", mf6_config_file)

# Objeto para acceder a toda la funcionalidad de la API
mf6 = ModflowApi(mf6_lib_path, working_directory=o_sim.sim_path)

# Inicialización del modelo
mf6.initialize(mf6_config_file)

# Para la solución obtenida con np.linalg.solve()
SOL = np.zeros(dis['ncol'])

# Obtenemos el tiempo actual y el tiempo final de la simulación
current_time = mf6.get_current_time()
end_time = mf6.get_end_time()

# Máximo número de iteraciones para el algorimo de solución numérica
max_iter = mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))

linea = 50*chr(0x2015)
print(linea)
print("Iniciando la simulación")
print(linea)

# Ciclo sobre tiempo
while current_time < end_time:
    # Obtenemos el paso de tiempo
    dt = mf6.get_time_step()
    print("dt:", dt, ", t:", current_time, ", end_t:", end_time, ", max_iter:", max_iter)

    # Preparar el objeto de la API para obtener la solución y
    # con el paso de tiempo
    mf6.prepare_time_step(dt)
    mf6.prepare_solve()
    
    # Ciclo del algoritmo numérico de solución
    kiter = 0
    while kiter < max_iter:
        print("\nkiter :", kiter)
        
        # Construye el sistema del problema y lo resuelve
        has_converged = mf6.solve(1)
        
        if has_converged:
            print(f" ---> ¿Convergencia obtenida? : {has_converged}")
            break
        else:
            print(f" ---> ¿Convergencia obtenida? : {has_converged}")
            
        kiter += 1

    # En este momento podemos construir la matriz del sistema
    A, _, _, _ = xmf6.api.build_mat(mf6)
    RHS = mf6.get_value(mf6.get_var_address("RHS", 'SLN_1'))
    print("\nA:\n", A)
    print("\nRHS:\n", RHS)

    # Calculamos la solución con np.linalg.solve() para comparar
    SOL[:] = np.linalg.solve(A, RHS)
    print("\nSOL:\n", SOL)
        
    # Finalizamos la solución del paso de tiempo actual
    mf6.finalize_solve()

    # Finalizamos el paso de tiempo actual. 
    mf6.finalize_time_step()

    # Avanzamos en el tiempo
    current_time = mf6.get_current_time()

    if not has_converged:
        print("model did not converge")
        break

# Almacenamos la solución obtenida por MF6 (ojo: necesitamos hacer una copia del arreglo)
SOL_MF6 = np.copy(mf6.get_value_ptr(mf6.get_var_address("X", 'FLOW')))

# Finalizamos la simulación completa
try:
    mf6.finalize()
    success = True
except:
    raise RuntimeError

print(linea)
print("SOl (MF6):", SOL_MF6)
print(linea)
print("Finalizando la simulación")
print(linea)

###########


# --- Recuperamos los resultados de la simulación ---
head = xmf6.gwf.get_head(o_gwf)

x, y, z = o_gwf.modelgrid.xyzcellcenters
row_length = o_gwf.modelgrid.extent[1]

# --- Definición de la figura ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize =(8,4), height_ratios=[2,0.1])

# --- Gráfica 1. ---
iskip = 3
# Carga hidráulica obtenida con mf6.solve() almacenada en memoria.
ax1.scatter(x[0][::iskip], SOL_MF6[::iskip], marker="o", s = 40, fc = 'darkorange', ec = 'k', alpha=0.75, label="Mf6", zorder=5)

# Carga hidráulica obtenida con np.linalg.solve() almacenada en memoria.
ax1.scatter(x[0][::iskip], SOL[::iskip], marker="^", s = 10, c='blue', alpha = 0.75, label='np.linalg.solve()', zorder=5)

ax1.plot(x[0], SOL, c="k", zorder=0)

ax1.set_xlim(0, 12)
ax1.set_xticks(ticks=np.linspace(0, row_length,13))
ax1.set_xlabel("Distance (cm)")
ax1.set_ylabel("Head")
ax1.grid()
ax1.legend()

# --- Gráfica 2. ---
pmv = flopy.plot.PlotMapView(o_gwf, ax=ax2)
pmv.plot_grid(colors='dimgray', lw=0.5)
ax2.set_yticks(ticks=[0, 0.1])#, fontsize=8)
ax2.set_title("Mesh")

plt.tight_layout()
plt.show()