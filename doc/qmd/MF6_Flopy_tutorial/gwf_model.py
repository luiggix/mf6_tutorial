"""
@author: Luis M. de la Cruz [Updated on Wed May 14 18:00:33 CST 2025].
"""

import flopy

def build(conf, time, mesh, ic_data, chd_data, k_data):
    """
    Genera un modelo GWF con malla, condiciones de frontera y permeabilidad hidráulica adaptables.

    Parameters
    ----------
    conf: dict
        Diccionario con el nombre de la simulación, el nombre del ejecutable, la ruta del espacio de trabajo.

    time: dict
        Diccionario con los datos para la discretización temporal.

    mesh: dict
        Diccionario con los datos de la discretización espacial.

    ic_data: dict
        Diccionario con los datos de las condiciones iniciales.
        
    chd_data: numpy.ndarray
        Datos de las condiciones de frontera en una malla rectangular.

    k_data: numpy.ndarray
        Datos de la permeabilidad hidráulica.
    
    Returns
    -------
    o_sim:
        Objeto de la simulación.
        
    o_gwf:
        Objeto del modelod GWF.
    """
    # --- Inicialización de la simulación ---
    o_sim = flopy.mf6.MFSimulation(
        sim_name = conf["sim_name"],
        exe_name = conf["exe_name"],
        sim_ws   = conf["sim_ws"]
    )

    ### --- COMPONENTES ---

    # --- Tiempos de simulación ---
    o_tdis = flopy.mf6.ModflowTdis(
        simulation = o_sim,
        time_units = time["units"],
        nper = time["nper"],
        perioddata = time["perioddata"]
    )
    
    # --- Solucionador numérico (IMS) ---
    o_ims = flopy.mf6.ModflowIms(
        simulation = o_sim,
    )
    
    # --- Modelo GWF ---
    o_gwf = flopy.mf6.ModflowGwf(
        simulation = o_sim,
        modelname = conf["sim_name"],
        model_nam_file = f"{conf["sim_name"]}.nam",
        save_flows = True # Almacena los flujos, particularmente el budget
    )
    
    ### --- PAQUETES ---

    # --- Discretización espacial ---º
    o_dis = flopy.mf6.ModflowGwfdis(
        model = o_gwf,
        length_units = "METERS",
        nlay = mesh['nlay'],
        nrow = mesh['nrow'],
        ncol = mesh['ncol'],
        delr = mesh['delr'],
        delc = mesh['delc'],
        top  = mesh['top'],
        botm = mesh['botm'],
    )
    
    # --- Condiciones iniciales ---
    o_ic = flopy.mf6.ModflowGwfic(
        model = o_gwf,
        strt = ic_data["strt"]
    )
    
    # --- Condiciones de frontera ---
    o_chd = flopy.mf6.ModflowGwfchd(
        model = o_gwf,
        stress_period_data=chd_data,
    )
    
    # --- Propiedades de flujo ---
    o_npf = flopy.mf6.ModflowGwfnpf(
        model = o_gwf,
        save_specific_discharge = True,
        k = k_data, # conductividad hidráulica
    )
    
    # --- Configuración de la salida ---
    o_oc = flopy.mf6.ModflowGwfoc(
        model = o_gwf,
        budget_filerecord = f"{conf["sim_name"]}.bud",
        head_filerecord = f"{conf["sim_name"]}.hds",
        saverecord = [("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord = [("HEAD", "ALL")]
    )

    return o_sim, o_gwf