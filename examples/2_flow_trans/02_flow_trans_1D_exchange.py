import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
import xmf6

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
    
phys = dict(
    specific_discharge = 0.1,  # Specific discharge ($cm s^{-1}$)
    hydraulic_conductivity = 0.01,  # Hydraulic conductivity ($cm s^{-1}$)
    source_concentration = 1.0,  # Source concentration (unitless)
    porosity = 0.1,  # Porosity of mobile domain (unitless)
    initial_concentration = 0.0,  # Initial concentration (unitless)
)
xmf6.nice_print(phys, "Parámetros físicos")


# Parámetros de la simulación (flopy.mf6.MFSimulation)
init = {
    'sim_name' : "flow",
#    'exe_name' : "C:\\Users\\luiggi\\Documents\\GitSites\\mf6_tutorial\\mf6\\windows\\mf6",
    'exe_name' : "../../mf6/macosarm/mf6",
    'sim_ws' : "output_ft_ex"
}

# Parámetros para el tiempo (flopy.mf6.ModflowTdis)
tdis = {
    'units': "seconds",
    'nper' : 1,
    'perioddata': [(120.0, 240, 1.0)]
}

# Parámetros para la solución numérica (flopy.mf6.ModflowIms)
#ims = {}


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
o_sim = xmf6.common.init_sim(silent = True, init = init, tdis = tdis)
o_gwf, packages = xmf6.gwf.set_packages(o_sim, silent = True,
                                        gwf = gwf, 
                                        dis = dis, ic = ic, chd = chd, npf = npf, oc = oc, well = well)

o_ims = flopy.mf6.ModflowIms(
    o_sim,
    print_option="ALL",
#        outer_dvclose=hclose,
#        outer_maximum=nouter,
    under_relaxation="NONE",
#        inner_maximum=ninner,
#        inner_dvclose=hclose,
#        rcloserecord=rclose,
    linear_acceleration="BICGSTAB",
    scaling_method="NONE",
    reordering_method="NONE",
#        relaxation_factor=relax,
    filename=f"{o_gwf.name}.ims",
)
o_sim.register_ims_package(o_ims, [o_gwf.name])

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
    'sim_name' : "trans",
#    'exe_name' : "C:\\Users\\luiggi\\Documents\\GitSites\\mf6_tutorial\\mf6\\windows\\mf6",
    'exe_name' : "../../mf6/macosarm/mf6",
##    'sim_ws' : "output_flow_trans_1D"
}

# Parámetros para el tiempo (flopy.mf6.ModflowTdis)
#El mismo tdis de flow

# Parámetros para la solución numérica (flopy.mf6.ModflowIms)
#ims_t = {
#    'linear_acceleration': "bicgstab"
#}

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
#fmi = {
#    "packagedata" : [("GWFHEAD", f"{init['sim_name']}.hds", None),
#                     ("GWFBUDGET", f"{init['sim_name']}.bud", None),
#                    ]
#}

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
    'budget_filerecord': f"{init_t['sim_name']}.cbc",
    'concentration_filerecord': f"{init_t['sim_name']}.ucn",
    'saverecord' : [("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
    'printrecord' : [("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
}

# --- Inicialización de la simulación ---
#o_sim_t = xmf6.common.init_sim(silent = True, init = init_t, 
#                               tdis = tdis, ims = ims_t)
o_gwt, packagest = xmf6.gwt.set_packages(o_sim, silent = True,
                                        gwt = gwt, 
                                        dis = dis, ic = ic_t, 
                                        mst = mst, adv = adv, dsp = dsp, ssm = ssm, 
                                        oc = oc_t)

imsgwt = flopy.mf6.ModflowIms(
        o_sim,
        print_option="ALL",
#        outer_dvclose=hclose,
#        outer_maximum=nouter,
        under_relaxation="NONE",
#        inner_maximum=ninner,
#        inner_dvclose=hclose,
#        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
#        relaxation_factor=relax,
        filename=f"{o_gwt.name}.ims",
)
o_sim.register_ims_package(imsgwt, [o_gwt.name])

o_obs = xmf6.common.set_obs(o_gwt, obs, silent = True)

### Aquí va el intercambio
flopy.mf6.ModflowGwfgwt(
    o_sim, exgtype="GWF6-GWT6", exgmnamea=o_gwf.name, exgmnameb=o_gwt.name
)

# --- Escritura de archivos ---
o_sim.write_simulation(silent = True)

# --- Ejecución de la simulación ---
o_sim.run_simulation(silent = False)


# --- Recuperamos los resultados de flujo de la simulación ---
head = xmf6.gwf.get_head(o_gwf)
qx, qy, qz, n_q = xmf6.gwf.get_specific_discharge(o_gwf, text="DATA-SPDIS")


# --- Recuperamos los resultados de la concentración de la simulación ---
cobj = o_gwt.output.concentration()
times = cobj.get_times()
times = np.array(times)

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
#    gwt_conc = xmf6.gwt.get_concentration(o_sim, t)
    conc = cobj.get_data(totim=t)
    ax3.scatter(x[0][::iskip], conc[0, 0][::iskip], label=f'GWT. Time = {t}',
                ec="k", alpha=0.75, s=20, zorder=5)    
    
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
