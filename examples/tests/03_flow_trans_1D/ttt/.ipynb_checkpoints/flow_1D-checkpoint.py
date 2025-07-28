import flopy

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


if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt # Graficación
    import os
    
    from mesh import MeshDis
    from tdis import TDis
    from vis import plot1D
    from output import nice_print, OFiles
    
    # ----- Definición de Parámetros -----
    mesh = MeshDis(
        nrow = 1,    # Number of rows
        ncol = 120,  # Number of columns
        nlay = 1,    # Number of layers
        row_length = 12.0,    # Length of system ($cm$)
        column_length = 0.1,  # Length of system ($cm$)
        top = 1.0,   # Top of the model ($cm$)
        bottom = 0,  # Layer bottom elevation ($cm$)
    )
    nice_print(mesh.get_dict(), 'Space discretization')

    tdis = TDis(
        perioddata = ((120, 240, 1.0),) # PERLEN, NSTP, TSMULT
    )
    nice_print(tdis, 'Time discretization') 
    
    ml_units = {
        "time": "seconds",
        "length": "centimeters"
    }
    nice_print(ml_units, 'Units')

    ph_par = dict(
        specific_discharge = 0.1,  # Specific discharge ($cm s^{-1}$)
        hydraulic_conductivity = 0.01,  # Hydraulic conductivity ($cm s^{-1}$)
        source_concentration = 1.0  # Source concentration (unitless)
    )
    nice_print(ph_par, 'Physical parameters')

    os_par = dict(
        ws = os.getcwd() + '/output', # Ruta de donde estamos actualmente
        mf6_exe = '/home/jovyan/GMMC/WMA/mf6/bin/mf6', # Ejecutable
        flow_name = 'flow', # Nombre de la simulación             
    )
    nice_print(os_par, 'MODFLOW 6 environment')
                     
    oc_par = dict(
        head_file = f"{os_par['flow_name']}.hds", 
        fbudget_file = f"{os_par['flow_name']}.bud",            
    )
    nice_print(oc_par, 'Output files')
    # ------------------------------------
    
    sim, gwf = build_gwf_1D(mesh, tdis, ph_par, ml_units, os_par, oc_par)
    sim.write_simulation()
    sim.run_simulation()

    of = OFiles(gwf, os_par, oc_par)
    head = of.get_head()

    xi, _, _ = mesh.get_coords()
    plt.figure(figsize=(10,3))
    ax = plt.gca()
    ax.set_xlim(0, 12)
    ax.set_xticks(ticks=np.linspace(0, mesh.row_length,13))
    ax.set_xlabel("Distance (cm)")
    ax.set_ylabel("Head (unitless)")
    ax.grid(True)
    plt.tight_layout()
    plot1D(ax, xi, head[0, 0], marker=".", ls ="-", mec="blue", mfc="none", markersize="1", label = 'Head')
