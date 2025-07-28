
import os, sys   # Intefaces con el sistema operativo.
import numpy as np # Manejo de arreglos numéricos multidimensionales
import matplotlib.pyplot as plt # Graficación

# Biblioteca y módulos de flopy
import flopy
from flopy.plot.styles import styles

# Definición de parámetros de decaimiento y sorción
from sorption_decay import *
from flow_1D import build_gwf_1D
import xmf6

def build_gwt_1D(mesh, tdis, ph_par, ml_units, os_par, oc_par):
    """
    Función que crea una simulación para transport usando GWT

    Parameters
    ----------
    mesh: MeshDis
    Objeto que gestiona atributos y métodos de una malla rectangular
    estructurada y uniforme.

    tdis: TDis
    Objeto que gestiona atributos del tiempo.

    ph_par: dict
    Diccionario con los parámetros físicos del problema.

    ml_units: dict
    Unidades para los parámetros del problema.

    os_par: dict
    Parámetros para ejecución de MODFLOW 6.
    
    oc_par: dict
    Parámetros para los archivos de salida.

    Returns
    -------
    sim: MFSimulation
    Objeto de la simulación.

    gwf:
    Objeto del modelo GWF.
    """
    # Creamos la simulación
    sim = flopy.mf6.MFSimulation(
        sim_name=os_par["tran_name"], 
        sim_ws=os_par["ws"], 
        exe_name=os_par["mf6_exe"]
    )

    # Definimos la componente para el tiempo
    flopy.mf6.ModflowTdis(
        sim, 
        nper=tdis.nper(), 
        perioddata=tdis.perioddata(), 
        time_units=ml_units["time"]
    )

    # Definimos la componente para la solución numérica
    flopy.mf6.ModflowIms(
        sim,
        linear_acceleration="bicgstab"        
    )

    # Definimos el modelo GWT
    gwt = flopy.mf6.ModflowGwt(
        sim,
        modelname=os_par["tran_name"],
        save_flows=True
    )
    
    # Paquete para discretización espacial
    flopy.mf6.ModflowGwfdis(
        gwt,
        length_units=ml_units["length"],
        nlay=mesh.nlay,
        nrow=mesh.nrow,
        ncol=mesh.ncol,
        delr=mesh.delr,
        delc=mesh.delc,
        top=mesh.top,
        botm=mesh.bottom,
    )

    # Paquete para las condiciones iniciales
    flopy.mf6.ModflowGwfic(
        gwt, 
        strt=0.0
    )

    # Paquete para ...
    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=ph_par["porosity"],
        **get_sorption_dict(ph_par["retardation_factor"]),
        **get_decay_dict(ph_par["decay_rate"], 
                         ph_par["retardation_factor"] > 1.0),
    )

    # Paquete para la adveccción
    flopy.mf6.ModflowGwtadv(
        gwt, 
        scheme="TVD"
    )

    # Paquete para ...
    flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=True,
        alh=ph_par["longitudinal_dispersivity"],
        ath1=ph_par["longitudinal_dispersivity"],
    )

    # Paquete para ...
    pd = [ 
        ("GWFHEAD", oc_par["head_file"], None),
        ("GWFBUDGET", oc_par["fbudget_file"], None),
    ] 
    flopy.mf6.ModflowGwtfmi(
        gwt, 
        packagedata=pd
    )

    # Paquete para ...
    sourcerecarray = [["WEL-1", "AUX", "CONCENTRATION"]]
    flopy.mf6.ModflowGwtssm(
        gwt, 
        sources=sourcerecarray
    )

    # Paquete para ...
    obs_data = {
        "transporte.obs.csv": [
            ("X005", "CONCENTRATION", (0, 0, 0)),
            ("X405", "CONCENTRATION", (0, 0, 40)),
            ("X1105", "CONCENTRATION", (0, 0, 110)),
        ],
    }
    
    flopy.mf6.ModflowUtlobs(
        gwt, 
        digits=10, 
        print_input=True, 
        continuous=obs_data
    )

    # Paquete para la salida
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=oc_par["tbudget_file"],
        concentration_filerecord=oc_par["concentration_file"],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, gwt

from wexler1_test import sol_analytical_t
def plot_tran_1D(sim, mesh, tdis, ph_par, os_par, oc_par, savefig = False):
    """
    Función para graficar los resultados.

    Paramaters
    ----------
    gwf: ModflowGwf
    Objeto del modelo de flujo GWF

    mesh: MeshDis
    Objeto que gestiona atributos y métodos de una malla rectangular
    estructurada y uniforme.

    os_par: dict
    Parámetros para ejecución de MODFLOW 6, archivos de salida y 
    path del workspace.
    """
    mf6gwt_ra = sim.get_model("transport").obs.output.obs().data
    ucnobj_mf6 = sim.transport.output.concentration()
    simtimes = mf6gwt_ra["totim"]
    obsnames = ["X005", "X405", "X1105"]

    with styles.USGSPlot():
        plt.rcParams['font.family'] = 'DeJavu Sans'
    
        fig, axs = plt.subplots(2, 1, figsize=(5,6), tight_layout=True)
    
        iskip = 5
    
        atimes = np.arange(0, tdis.total_time(), 0.1)
        
        for i, x in enumerate([0.05, 4.05, 11.05]):
            a1, idx_filter = sol_analytical_t(i, x, atimes,mesh, ph_par) 
                    
            axs[0].plot(atimes[idx_filter], a1[idx_filter], color="k", label="ANALYTICAL")
    
            axs[0].plot(simtimes[::iskip], mf6gwt_ra[obsnames[i]][::iskip],
                        marker="o", ls="none", mec="blue", mfc="none", markersize="4",
                        label="MODFLOW 6")
            axs[0].set_ylim(-0.05, 1.2)
            axs[0].set_xlim(0, 120)
            axs[0].set_xlabel("Time (seconds)")
            axs[0].set_ylabel("Normalized Concentration (unitless)")
            
        ctimes = [6.0, 60.0, 120.0]
        x, _, _ = mesh.get_coords()
        for i, t in enumerate(ctimes):
            a1, idx_filter = sol_analytical_t(i, x, t, mesh, ph_par, False)
            
            axs[1].plot(x, a1, color="k", label="ANALYTICAL")
            simconc = ucnobj_mf6.get_data(totim=t).flatten()
            axs[1].plot(x[::iskip], simconc[::iskip],
                        marker="o", ls="none", mec="blue", mfc="none", markersize="4",
                        label="MODFLOW 6")
            axs[1].set_ylim(0, 1.1)
            axs[1].set_xlim(0, 12)
            axs[1].set_xlabel("Distance (cm)")
            axs[1].set_ylabel("Normalized Concentration (unitless)")
        
        if savefig:
            plt.savefig('conc.pdf')
        else:
            plt.show()


if __name__ == '__main__':

    # ----- Definición de Parámetros -----
    mesh = xmf6.MeshDis(
        nrow = 1,    # Number of rows
        ncol = 120,  # Number of columns
        nlay = 1,    # Number of layers
        row_length = 12.0,    # Length of system ($cm$)
        column_length = 0.1,  # Length of system ($cm$)
        top = 1.0,   # Top of the model ($cm$)
        bottom = 0,  # Layer bottom elevation ($cm$)
    )
    xmf6.nice_print(mesh.get_dict(), 'Space discretization')

    tm_par = dict(
        nper = 1,  # Number of periods
        total_time = 120.0,  # Simulation time ($s$)
        nstp = 240,   # Number of time steps
        tsmult = 1.0  # Multiplier for the length of successive time steps.
    )
    xmf6.nice_print(tm_par, 'Time discretization')
    
    ml_units = {
        "time": "seconds",
        "length": "centimeters"
    }
    xmf6.nice_print(ml_units, 'Units')


    ph_par = dict(
        specific_discharge = 0.1,  # Specific discharge ($cm s^{-1}$)
        hydraulic_conductivity = 0.01,  # Hydraulic conductivity ($cm s^{-1}$)
        source_concentration = 1.0,  # Source concentration (unitless)
        porosity = 0.1,  # Porosity of mobile domain (unitless)
        initial_concentration = 0.0,  # Initial concentration (unitless)
        longitudinal_dispersivity = 0.1,
        retardation_factor = 1.0,
        decay_rate =  0.0
    )
    ph_par["dispersion_coefficient"] = ph_par["longitudinal_dispersivity"] * ph_par["specific_discharge"] / ph_par["retardation_factor"]
    xmf6.nice_print(ph_par, 'Physical parameters')
    
    os_par = dict(
        ws = os.getcwd(), # Ruta de donde estamos actualmente
        mf6_exe = '/home/jovyan/GMMC/WMA/mf6/bin/mf6', # Ejecutable
        flow_name = 'flow', # Nombre de la simulación para flujo
        tran_name = 'transport' # Nombre de la simulación para transporte
    )
    xmf6.nice_print(os_par, 'MODFLOW 6 environment')
                     
    oc_par = dict(
        head_file = f"{os_par['flow_name']}.hds",
        fbudget_file = f"{os_par['flow_name']}.bud",
        concentration_file=f"{os_par['tran_name']}.ucn",
        tbudget_file = f"{os_par['tran_name']}.bud",
    )
    xmf6.nice_print(oc_par, 'Output files')
    # ------------------------------------

    # Solución de flujo
    sim_f, gwf = build_gwf_1D(mesh, tm_par, ph_par, ml_units, os_par, oc_par)
    sim_f.write_simulation(silent=True)
    sim_f.run_simulation(silent=True)
    xmf6.plot_flow_1D(gwf, mesh, os_par, oc_par, True)

    # Solución de transporte
    sim_t, gwt = build_gwt_1D(mesh, tm_par, ph_par, ml_units, os_par, oc_par)
    sim_t.write_simulation()
    sim_t.run_simulation()
    xmf6.plot_tran_1D(sim_t, mesh, tm_par, ph_par, os_par, oc_par, True)
    