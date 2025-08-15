import numpy as np
import matplotlib.pyplot as plt
import os, sys
import flopy
from modflowapi import ModflowApi
import xmf6

def build_mat(mf6):
    """
    Construye la matriz del sistema.

    Parameters
    ----------
    mf6: ModflowApi
        Objeto para accedar a toda la funcionalidad de la API.
    """
    # Obtiene el número de renglones y columnas del sistema
    NCOL = mf6.get_value(mf6.get_var_address("NCOL", 'SLN_1'))
    NROW = mf6.get_value(mf6.get_var_address("NROW", 'SLN_1'))

    # Obtiene los coeficientes de la matriz en formato CRS (Compressed Row Storage)
    # A: Coeficientes, JA: índices de la columna, IA: índice de inicio del renglón en JA.
    A = mf6.get_value(mf6.get_var_address("AMAT", 'SLN_1'))
    IA = mf6.get_value(mf6.get_var_address("IA", 'SLN_1'))
    JA = mf6.get_value(mf6.get_var_address("JA", 'SLN_1'))

    # Arreglo para almacenar la matriz en formato completo.
    Atest = np.zeros((NROW[0], NCOL[0]))
    idx = 0
    i = 0
    istart = IA[0] # Inicio del renglón en IA
    for iend in IA[1:]: # Recorremos desde el inicio de cada renglón
        for j in range(istart, iend): # Recorremos todos los elementos del renglón
            Atest[idx, JA[j-1]-1] = A[i] # Agregamos el coeficiente en la matriz completa
            i += 1
        istart = iend
        idx += 1
    return Atest, A, IA, JA # Regresamos la matriz densa y en el format CRS

# --- Preparación de la simulación ---

# --- Componentes ---

# Parámetros de la simulación (flopy.mf6.MFSimulation)
init = {
    'sim_name' : "flow",
    'exe_name' : "C:\\Users\\luiggi\\Documents\\GitSites\\mf6_tutorial\\mf6\\windows\\mf6",
#    'exe_name' : "../../mf6/macosarm/mf6",
    'sim_ws' : "sandbox4"
}

# Parámetros para el tiempo (flopy.mf6.ModflowTdis)
tdis = {
    'units': "DAYS",
    'nper' : 3,
    'perioddata': [(1.0, 5, 1.0), (2.0, 3, 1.2), (3.0, 4, 1.1)] 
}

# Parámetros para la solución numérica (flopy.mf6.ModflowIms)
ims = {}

# Parámetros para el modelo de flujo (flopy.mf6.ModflowGwf)
gwf = { 
    'modelname': init["sim_name"],
    'model_nam_file': f"{init["sim_name"]}.nam",
    'save_flows': True
}

# --- Paquetes del modelo de flujo ---
lx = 25
ly = 25
nrow = 20 #int(input("nrow = "))
ncol = 20 #int(input("ncol = "))
#if nrow > 10 or nrow < 0 or ncol > 10 or ncol < 0:
#    print("Este ejemplo está limitado a:")
#    print(" 0 < nrow =< 10")
#    print(" 0 < ncol =< 10")
#    sys.exit()
delr = lx / ncol
delc = ly / nrow 

# Parámetros para la discretización espacial (flopy.mf6.ModflowGwfdis)
dis = {
    'length_units': "meters",
    'nlay': 5, 
    'nrow': nrow, 
    'ncol': ncol,
    'delr': delr, 
    'delc': delc, 
    'top' : 0.0, 
    'botm': [-2.0, -4.0, -6.0, -8.0, -10.0] 
}

# Parámetros para las condiciones iniciales (flopy.mf6.ModflowGwfic)
ic = {
    'strt': 1.0
}

# Parámetros para las condiciones de frontera (flopy.mf6.ModflowGwfchd)
chd_data = []
for row in range(dis['nrow']):
    chd_data.append([(0, row, 0), 10.0])       # Condición en la pared izquierda
    chd_data.append([(0, row, dis['ncol'] - 1), 5.0]) # Condición en la pared derecha

chd = {
    'stress_period_data': chd_data,     
}

# Parámetros para las propiedades de flujo (flopy.mf6.ModflowGwfnpf)
npf = {
    'save_specific_discharge': True,
    'save_saturation' : True,
    'icelltype' : 0,
    'k' : 0.01,
}


# === Pozo en el centro con caudal variable ===
well_data = {
    0: [[(1, 5, 5), -0.5]],
    1: [[(1, 5, 5), -0.025]],
    2: [[(1, 5, 5), -0.015]],
}

well = {
    'stress_period_data': well_data,
    'pname': "WEL-1",
    'save_flows': True
}

# Parámetros para almacenar y mostrar la salida de la simulación (flopy.mf6.ModflowGwfoc)
oc = {
    'budget_filerecord': f"{init['sim_name']}.bud",
    'head_filerecord': f"{init['sim_name']}.hds",
    'saverecord': [("HEAD", "ALL"), ("BUDGET", "ALL")],
    'printrecord': [("HEAD", "ALL")]
}

# Inicialización de la simulación
o_sim = xmf6.common.init_sim(init = init, tdis = tdis, ims = ims, silent = True)

# Configuración de los paquetes para el modelo de flujo
o_gwf, package_list = xmf6.gwf.set_packages(o_sim, silent = True,
                                            gwf = gwf, dis = dis, ic = ic, chd = chd, npf = npf, well = well, oc = oc)

# Escritura de los archivos de entrada
o_sim.write_simulation(silent = False)

# Ejecución de la simulación
o_sim.run_simulation()

# --- Recuperamos los resultados de la simulación ---
head, harray = xmf6.gwf.get_head(o_gwf, binary=True)
qx, qy, qz, n_q = xmf6.gwf.get_specific_discharge(o_gwf, text="DATA-SPDIS")

print(harray.shape)
# Obtener tiempos disponibles
times = head.get_times()
print("Tiempos disponibles:")
[print(t) for t in times]

# Obtener número de periodos y pasos
nper = len(times)
print(nper)

# --- Parámetros para las gráficas ---
grid = o_gwf.modelgrid
x, y, z = grid.xyzcellcenters
hvmin = np.nanmin(harray)
hvmax = np.nanmax(harray)
xticks = x[0]
yticks = y[:,0]
zticks = z[:,0,0]
xlabels = [f'{x:1.1f}' for x in xticks]
ylabels = [f'{y:1.1f}' for y in yticks]
zlabels = [f'{z:1.1f}' for z in zticks]

ilay = int(input("Capa [1,2,3,4,5] = ")) - 1

# --- Definición de la figura ---
fig, (ax1, ax2) = plt.subplots(1, 2, width_ratios=[1,1])

# --- Gráfica 1. ---
hview = flopy.plot.PlotMapView(model = o_gwf, ax = ax1)
hview.plot_grid(linewidths = 0.5, alpha = 0.5)
h_ac = hview.plot_array(harray[ilay,:,:], cmap = "YlGnBu", vmin = hvmin, vmax = hvmax, alpha = 0.75)
h_cb = plt.colorbar(h_ac, ax = ax1, label = "", cax = xmf6.vis.cax(ax1, h_ac))
h_cb.ax.tick_params(labelsize=6)

# Visualización de los vectores de la descarga específica.
ax1.quiver(x, y, qx[0], qy[0], scale = 0.05, 
           color = 'k', linewidth = 0.95, pivot = 'middle')

ax1.set_title("$h$ (Mf6)", fontsize=10)
ax1.set_ylabel("$y$ (m)", fontsize = 8)
ax2.set_xlabel("$x$ (m)", fontsize = 8)
ax1.set_xticks(ticks = [], labels = [], fontsize=6, rotation=90)
ax1.set_yticks(ticks = yticks, labels = yticks, fontsize=6)
ax1.set_xticks(ticks = xticks, labels = xticks, fontsize=6, rotation=90)
ax1.set_aspect('equal')

# --- Gráfica 2. ---
xz6 = flopy.plot.PlotCrossSection(model=o_gwf, ax = ax2, line={"row": 3})
xz6.plot_grid(linewidths = 0.5, alpha = 0.5)
xz6_ac = xz6.plot_array(harray, cmap = "YlGnBu", vmin = hvmin, vmax = hvmax, alpha = 0.75)
xz6_cb = plt.colorbar(xz6_ac, ax = ax2, label = " ", cax = xmf6.vis.cax(ax2, xz6_ac))
xz6_cb.ax.tick_params(labelsize=6)
ax2.set_ylabel("$z$ (m)", fontsize = 8)
ax2.set_xlabel("$x$ (m)", fontsize = 8)
ax2.set_xticks(ticks = xticks, labels = xticks, fontsize=6, rotation=90)
ax2.set_yticks(ticks = zticks, labels = zticks, fontsize=6)
ax2.set_aspect('equal')

plt.suptitle(f"Capa : {ilay+1}")
plt.tight_layout()
plt.show()
