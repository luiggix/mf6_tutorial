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
    porosity = 0.1,  # Porosity of mobile domain (unitless)
    initial_concentration = 0.0,  # Initial concentration (unitless)
)
xmf6.nice_print(phys, "Parámetros físicos")

############# GWF #############

# Parámetros de la simulación (flopy.mf6.MFSimulation)
init = {
    'sim_name' : "flow",
    'exe_name' : "C:\\Users\\luiggi\\Documents\\GitSites\\mf6_tutorial\\mf6\\windows\\mf6",
#    'exe_name' : "../../mf6/macosarm/mf6",
    'sim_ws' : "sandbox5"
}

# Parámetros para el tiempo (flopy.mf6.ModflowTdis)
tdis = {
    'units': "seconds",
    'nper' : 1,
    'perioddata': [(120.0, 240, 1.0)]
}

# Parámetros para la solución numérica (flopy.mf6.ModflowIms)
ims = {}

# Parámetros para el modelo de flujo (flopy.mf6.ModflowGwf)
gwf = { 
    'modelname': init["sim_name"],
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
    'k': phys["hydraulic_conductivity"] # Hydraulic conductivity ($cm s^{-1}$) 
}

# Parámetros para las propiedades de los pozos (flopy.mf6.ModflowGwfwel)
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
o_simf = xmf6.common.init_sim(silent = True, init = init, tdis = tdis, ims = ims)
o_gwf, packages = xmf6.gwf.set_packages(o_simf, silent = True,
                                        gwf = gwf, 
                                        dis = dis, ic = ic, chd = chd, npf = npf, oc = oc, well = well)
# --- Escritura de archivos ---
o_simf.write_simulation(silent = True)

# --- Ejecución de la simulación ---
o_simf.run_simulation()


############# GWT #############

def get_mst_dict(phys):
    mst = dict(
        porosity = phys["porosity"],
        sorption = None,
        bulk_density = None,
        distcoef = None, 
        first_order_decay = None,
        decay = None,
        decay_sorbed = None
    )
    
    if phys["retardation_factor"] > 1.0:
        mst["sorption"] = "linear"
        mst["bulk_density"] = 1.0
        mst["distcoef"] = (phys["retardation_factor"] - 1.0) * phys["porosity"] / mst["bulk_density"]

    if phys["decay_rate"] != 0.0:
        mst["first_order_decay"] = True
        mst["decay"] = phys["decay_rate"]
        if phys["retardation_factor"] > 1.0:
            mst["decay_sorbed"] = phys["decay_rate"]

    return mst

long_disp = [0.1, 1.0, 1.0, 1.0]
reta_fact = [1.0, 1.0, 2.0, 1.0]
deca_rate = [0.0, 0.0, 0.0, 0.01]
dir_names = ['p01a','p01b','p01c','p01d']

case = int(input("Caso (1-4)= ")) - 1 # 0, 1, 2, 3
dirname = dir_names[case]

# Agregamos más parámetros físicos al diccionario phys
phys["longitudinal_dispersivity"] = long_disp[case]
phys["retardation_factor"] = reta_fact[case]
phys["decay_rate"] = deca_rate[case]
phys["dispersion_coefficient"] = phys["longitudinal_dispersivity"] * phys["specific_discharge"] / phys["retardation_factor"]

xmf6.nice_print(phys, "Parámetros físicos")
print("Caso: {}".format(dirname))


# Parámetros de la simulación (flopy.mf6.MFSimulation)
init_t = {
    'sim_name' : "transport",
    'exe_name' : "C:\\Users\\luiggi\\Documents\\GitSites\\mf6_tutorial\\mf6\\windows\\mf6",
#    'exe_name' : "../../mf6/macosarm/mf6",
    'sim_ws' : "sandbox5"
}

# Parámetros para el tiempo (flopy.mf6.ModflowTdis)
#El mismo tdis de flow

# Parámetros para la solución numérica (flopy.mf6.ModflowIms)
ims_t = {
    'linear_acceleration': "bicgstab"
}

# Parámetros para el modelo de transporte (flopy.mf6.ModflowGwt)
gwt = { 
    'modelname': init_t["sim_name"],
    'save_flows': True
}

# Parámetros para la discretización espacial (flopy.mf6.ModflowGwtdis)
# El mismo dis de flow

# Parámetros para las condiciones iniciales (flopy.mf6.ModflowGwtic)
ic_t = {
    'strt': 0.0
}

# Parámetros para MST (flopy.mf6.ModflowGwtmst)
mst = get_mst_dict(phys)

# Parámetros para ADV (flopy.mf6.ModflowGwtadv)
adv = {
    "scheme" : "TVD"
}

# Parámetros para DSP (flopy.mf6.ModflowGwtdsp)
dsp = {
    "xt3d_off" : True,
    "alh" : phys["longitudinal_dispersivity"],
    "ath1" : phys["longitudinal_dispersivity"],
}

# Parámetros para FMI (flopy.mf6.ModflowGwtfmi)
fmi = {
    "packagedata" : [("GWFHEAD", f"{init['sim_name']}.hds", None),
                     ("GWFBUDGET", f"{init['sim_name']}.bud", None),
                    ]
}

# Parámetros para SSM (flopy.mf6.ModflowGwtssm)
ssm = {
    "sources" : [["WEL-1", "AUX", "CONCENTRATION"]]
}

# Parámetros para OBS (flopy.mf6.ModflowGwtobs)
obs = {
    "digits" : 10, 
    "print_input" : True, 
    "continuous" : {
        "transporte.obs.csv": [
            ("X005", "CONCENTRATION", (0, 0, 0)),
            ("X405", "CONCENTRATION", (0, 0, 40)),
            ("X1105", "CONCENTRATION", (0, 0, 110)),
        ],
    }
}

# Parámetros para almacenar y mostrar la salida de la simulación (flopy.mf6.ModflowGwtoc)
oc_t = {
    'budget_filerecord': f"{init_t['sim_name']}.bud",
    'concentration_filerecord': f"{init_t['sim_name']}.ucn",
    'saverecord' : [("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
    'printrecord' : [("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
}

# --- Inicialización de la simulación ---
o_sim_t = xmf6.common.init_sim(silent = True, init = init_t, 
                               tdis = tdis, ims = ims_t)
o_gwt, packagest = xmf6.gwt.set_packages(o_sim_t, silent = True,
                                        gwt = gwt, 
                                        dis = dis, ic = ic_t, 
                                        mst = mst, adv = adv, dsp = dsp, fmi = fmi, ssm = ssm, 
                                        oc = oc_t)

o_obs = xmf6.common.set_obs(o_gwt, obs, silent = True)

# --- Escritura de archivos ---
o_sim_t.write_simulation(silent = True)

# --- Ejecución de la simulación ---
#o_sim_t.run_simulation(silent = True)


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
mf6_config_file = os.path.join(o_sim_t.sim_path, 'mfsim.nam')
print("Shared library:", mf6_lib_path)
print("Config file:", mf6_config_file)

# Objeto para acceder a toda la funcionalidad de la API
mf6 = ModflowApi(mf6_lib_path, working_directory=o_sim_t.sim_path)

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
#    print("\nA:\n", A)
#    print("\nRHS:\n", RHS)

    # Calculamos la solución con np.linalg.solve() para comparar
    SOL[:] = np.linalg.solve(A, RHS)
#    print("\nSOL:\n", SOL)
        
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
SOL_MF6 = np.copy(mf6.get_value_ptr(mf6.get_var_address("X", 'TRANSPORT')))

# Finalizamos la simulación completa
try:
    mf6.finalize()
    success = True
except:
    raise RuntimeError

print(linea)
print("SOl:", SOL)
print(linea)
print("SOl (MF6):", SOL_MF6)
print(linea)
print("Finalizando la simulación")
print(linea)

###########

# --- Recuperamos los resultados de flujo de la simulación ---
head = xmf6.gwf.get_head(o_gwf)
qx, qy, qz, n_q = xmf6.gwf.get_specific_discharge(o_gwf, text="DATA-SPDIS")

# --- Recuperamos las coordenadas del dominio
x, y, z = o_gwf.modelgrid.xyzcellcenters
row_length = o_gwf.modelgrid.extent[1]

#--- Leemos la solución analítica de la ec. de transporte ---
sol_path = 'analytic/' + dirname
a1_0 = np.load(sol_path + '/a1_x_0.npy')
a1_1 = np.load(sol_path + '/a1_x_1.npy')
a1_2 = np.load(sol_path + '/a1_x_2.npy')

# --- Definición de la figura ---
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize =(6,6), height_ratios=[0.1,1,3])

# --- Gráfica 1. Malla ---
pmv = flopy.plot.PlotMapView(o_gwf, ax=ax1)
pmv.plot_grid(colors='dimgray', lw=0.5)
ax1.set_yticks(ticks=[0, 0.1])#, fontsize=8)
ax1.set_title("Mesh")

# --- Gráfica 2. Carga hidráulica---
ax2.plot(x[0], head[0, 0], marker=".", ls ="-", mec="blue", mfc="none", markersize="1", label = 'Head')
ax2.set_xlim(0, 12)
ax2.set_xticks(ticks=np.linspace(0, row_length,13))
#ax2.set_xlabel("Distance (cm)")
ax2.set_ylabel("Head")
ax2.grid()

# --- Gráfica 3. Concentración ---
# Solución analítica
ax3.plot(x[0], a1_0, c = 'k', label='Analytic')
ax3.plot(x[0], a1_1, c = 'k')
ax3.plot(x[0], a1_2, c = 'k')

# Solución numérica
citer = [11, 119, 239]
ctimes = [6.0, 60.0, 120.0]
iskip = 3
for c, (i, t) in enumerate(zip(citer, ctimes)):
    gwt_conc = xmf6.gwt.get_concentration(o_sim_t, t)
    ax3.scatter(x[0][::iskip], gwt_conc[::iskip], label=f'GWT. Time = {t}',
                ec="k", alpha=0.75, s=20, zorder=5)    

ax3.scatter(x[0][::iskip], SOL[::iskip], label=f'np.linalg.solve()',
            marker = "1", c="k", alpha=0.95, s=60, zorder=5)

#ax3.scatter(x[0][::iskip], SOL_MF6[::iskip], label=f'MF6',
#            marker = "v", ec="k", alpha=0.75, s=10, zorder=5)

# Decoración de la gráfica
ax3.legend(fontsize=9)
ax3.set_xlim(0, 12)
ax3.set_ylim(-0.1, 1.5)
ax3.set_xticks(ticks=np.linspace(0, row_length,13))
ax3.set_xlabel("Distance (cm)")
ax3.set_ylabel("Concentration")
ax3.grid(True)
ax3.set_title(f"Case:{dirname}")
plt.tight_layout()
plt.show()
