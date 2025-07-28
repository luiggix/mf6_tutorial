
import os, sys   # Intefaces con el sistema operativo.
import numpy as np # Manejo de arreglos numéricos multidimensionales
import matplotlib.pyplot as plt # Graficación

# Biblioteca y módulos de flopy
import flopy
from flopy.plot.styles import styles
import xmf6

def build_gwf_1D(mesh, tdis, ph_par, ml_units, os_par, oc_par):
    """
    Función que crea una simulación para flujo usando GWF

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
        sim_name=os_par["flow_name"], 
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
        sim
    )

    # Definimos el modelo GWF
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=os_par["flow_name"],
        save_flows=True
    )
    
    # Paquete para discretización espacial
    flopy.mf6.ModflowGwfdis(
        gwf,
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
        gwf, 
        strt=1.0
    )

    # Paquete para las propiedades de flujo en los nodos
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_saturation=True,
        icelltype=0,
        k=ph_par["hydraulic_conductivity"],
    )

    # Paquete CHD
    flopy.mf6.ModflowGwfchd(
        gwf, 
        stress_period_data=[[(0, 0, mesh.ncol - 1), 1.0]]
    ) 

    # Paquete de pozos
    q = ph_par["specific_discharge"] * mesh.delc * mesh.delr * mesh.top
    aux = ph_par["source_concentration"]
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=[[(0, 0, 0), q, aux,]],
        pname="WEL-1",
        auxiliary=["CONCENTRATION"],
    )

    # Paquete para la salida
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=oc_par["head_file"],
        budget_filerecord=oc_par["fbudget_file"],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, gwf

def plot_flow_1D(gwf, mesh, os_par, oc_par, savefig = False):
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
    # Obtenemos los resultados de la carga hidráulica
    head = flopy.utils.HeadFile(
        os.path.join(os_par["ws"], 
                     oc_par["head_file"])).get_data()

    print('head_L = {} \t head_R = {}'.format(head[0,0,0], head[0,0,-1]))
    
    # Obtenemos los resultados del BUDGET
    bud  = flopy.utils.CellBudgetFile(
        os.path.join(os_par["ws"], 
                     oc_par["fbudget_file"]),
        precision='double'
    )
    # Obtenemos las velocidades
    spdis = bud.get_data(text='DATA-SPDIS')[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)
    
    with styles.USGSPlot():
        plt.rcParams['font.family'] = 'DeJavu Sans'
        x, _, _ = mesh.get_coords()
        plt.figure(figsize=(10,3))
        plt.plot(x, head[0, 0], marker=".", ls ="-", mec="blue", mfc="none", markersize="1", label = 'Head')
        plt.xlim(0, 12)
        plt.xticks(ticks=np.linspace(0, mesh.row_length,13))
        plt.xlabel("Distance (cm)")
        plt.ylabel("Head (unitless)")
        plt.legend()
        plt.grid()

        if savefig:
            plt.savefig('head.pdf')
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
        nstp = 1.0,   # Number of time steps
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
        source_concentration = 1.0  # Source concentration (unitless)
    )
    xmf6.nice_print(ph_par, 'Physical parameters')
    
    os_par = dict(
        ws = os.getcwd(), # Ruta de donde estamos actualmente
        mf6_exe = '/home/jovyan/GMMC/WMA/mf6/bin/mf6', # Ejecutable
        flow_name = 'flow', # Nombre de la simulación
        head_file = "flow.hds", 
        fbudget_file = "flow.bud"
    )

    os_par = dict(
        ws = os.getcwd(), # Ruta de donde estamos actualmente
        mf6_exe = '/home/jovyan/GMMC/WMA/mf6/bin/mf6', # Ejecutable
        flow_name = 'flow', # Nombre de la simulación             
    )
    xmf6.nice_print(os_par, 'MODFLOW 6 environment')
                     
    oc_par = dict(
        head_file = f"{os_par['flow_name']}.hds", 
        fbudget_file = f"{os_par['flow_name']}.bud",            
    )
    xmf6.nice_print(oc_par, 'Output files')
    # ------------------------------------
    
    sim, gwf = build_gwf_1D(mesh, tm_par, ph_par, ml_units, os_par, oc_par)
    sim.write_simulation()
    sim.run_simulation()
    xmf6.plot_flow_1D(gwf, mesh, os_par, oc_par, True)


