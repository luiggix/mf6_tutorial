# exchange function between PRMS6 surface and soil BMIs
surf2soil_vars = ['hru_ppt', 'hru_area_perv', 'hru_frac_perv',
                  'infil', 'sroff', 'potet', 'hru_intcpevap',
                  'snow_evap', 'snowcov_area', 'soil_rechr', 'soil_rechr_max',
                  'soil_moist', 'soil_moist_max', 'hru_impervevap',
                  'srunoff_updated_soil', 'transp_on']

surf2soil_cond_vars = ['dprst_evap_hru', 'dprst_seep_hru']

soil_cond_vars = ['soil_rechr_chg', 'soil_moist_chg']

soil2surf_vars = ['infil', 'sroff', 'soil_rechr', 'soil_moist']

def soilinput(msurf, msoil, exch_vars, surf_cond_vars, soil_cond_vars,
              dprst_flag, dyn_dprst_flag, imperv_flag):
    for var in exch_vars:
        msoil.set_value(var, msurf.get_value(var))
    if dprst_flag == 1:
        for var in surf_cond_vars:
            msoil.set_value(var, msurf.get_value(var))
    if dyn_dprst_flag in [1, 3] or imperv_flag in [1, 3]:
        for var in soil_cond_vars:
            msoil.set_value(var, msurf.get_value(var))

def soil2surface(msoil, msurf, exch_vars):
    for var in exch_vars:
        msurf.set_value(var, msoil.get_value(var))

def update_coupled(msurf, msoil):
    dprst_flag = msurf.get_value('dprst_flag')
    dyn_dprst_flag = msoil.get_value('dyn_dprst_flag')
    dyn_imperv_flag = msoil.get_value('dyn_imperv_flag')
    msurf.update()
    soilinput(msurf, msoil, surf2soil_vars, surf2soil_cond_vars, soil_cond_vars,
              dprst_flag, dyn_dprst_flag, dyn_imperv_flag)
    msoil.update()
    soil2surface(msoil, msurf, soil2surf_vars)
