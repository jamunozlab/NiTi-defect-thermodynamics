#!/usr/bin/env python

"""
Extracts AMReX simulation data and puts it in ml-friendly np array files
@Author: Jorge Munoz
@Date: June 30, 2023; mods July 11th, July14th; additions July 17th; minor additions August 8th
"""

import numpy as np
import pandas as pd
import os, yt
yt.set_log_level(0)

def get_voltage_settings(simulation_path, input_filename='inputc'):
    with open(simulation_path+input_filename) as input_file:
        for line in input_file:
            split = line.split()
            you_get = split[0]
            if you_get == "Phi_Bc_hi":
                initial_voltage = float(split[-1])
            if you_get == "Phi_Bc_hi_max":
                final_voltage = float(split[-1])
            if you_get == "Phi_Bc_inc":
                voltage_step_increment = float(split[-1])
    
    #final_voltage is not guaranteed to be achieved in the simulation
    return initial_voltage, final_voltage, voltage_step_increment

def get_domain_settings(simulation_path, input_filename='inputc'):
    xyz_lo = []
    xyz_hi = []
    xyz_n_cell = []
    with open(simulation_path+input_filename) as input_file:
        for line in input_file:
            if "domain.prob_lo" in line:
                xyz_lo = [float(x) for x in line.split()[-3:]]
            if "domain.prob_hi" in line:
                xyz_hi = [float(x) for x in line.split()[-3:]]
            if "domain.n_cell" in line:
                xyz_n_cell = [int(x) for x in line.split()[-3:]]
                
    return xyz_lo, xyz_hi, xyz_n_cell

def get_geometry_settings(simulation_path, component='FE', input_filename='inputc'):
    lo = []
    hi = []
    with open(simulation_path+input_filename) as input_file:
        for line in input_file:
            if component + '_lo' in line:
                lo = [float(x) for x in line.split()[-3:]]
            if component + '_hi' in line:
                hi = [float(x) for x in line.split()[-3:]]
                
    return lo, hi

def get_scalar_from_input_file(scalar, simulation_path, input_filename='inputc'):
    with open(simulation_path+input_filename) as output_file:
        for line in output_file:
            if str(scalar) in line:
                return float(line.split()[-1])
        return None

def extract_voltage_changes_original(simulation_path, initial_voltage, voltage_step_increment, output_filename='output.txt'):
    #voltage = initial_voltage + voltage_step_increment
    voltage = initial_voltage
    voltage_list = [voltage]
    voltage_time_list = [0]
    with open(simulation_path+output_filename) as output_file:
        for line in output_file:
            if "voltage" in line:
                voltage_time_step = int(line.split()[-1])
                voltage_time_list.append(voltage_time_step)
                voltage = voltage + voltage_step_increment
                voltage_list.append(voltage)

    #print(pd.Series(voltage_list, index=voltage_time_list))
    return pd.Series(voltage_list, index=voltage_time_list)

def extract_voltage_changes(simulation_path, initial_voltage, voltage_step_increment, output_filename='output.txt'):
    #voltage = initial_voltage + voltage_step_increment
    
    voltage_list = []
    voltage_time_list = []
    
    for entry in os.scandir(simulation_path):
        if entry.is_dir() and 'plt' in entry.name:
            if entry.name == 'plt00000000': continue
            plt_name = entry.name
            voltage_time_step = int(plt_name[3:])
            voltage_time_list.append(voltage_time_step)
            #voltage = voltage + voltage_step_increment
            #voltage_list.append(voltage)
            
    voltage_time_series = pd.Series(voltage_time_list)
    voltage_time_series = voltage_time_series.sort_values()
    idx=voltage_time_series.values
    voltage_list = []
    voltage_time_list = []
    voltage = initial_voltage
    for i, ts in enumerate(voltage_time_series.values):
        voltage_list.append(voltage)
        voltage_time_list.append(ts)
        voltage = voltage + voltage_step_increment

    return pd.Series(voltage_list, index=voltage_time_list)

def etl(simulation_path, input_filename='inputc', output_filename='output.txt', verbose=False):
    initial_voltage, final_voltage, voltage_step_increment = get_voltage_settings(simulation_path, input_filename)
    #print(initial_voltage, final_voltage, voltage_step_increment)
    voltage_changes_series = extract_voltage_changes(simulation_path, initial_voltage=initial_voltage, voltage_step_increment=voltage_step_increment)
    
    dt = get_scalar_from_input_file('dt', simulation_path)
    
    xyz_lo, xyz_hi, xyz_n_cell = get_domain_settings(simulation_path, input_filename)
    
    epsilon_0 = get_scalar_from_input_file('epsilon_0', simulation_path)
    epsilon_DE = get_scalar_from_input_file('epsilon_de', simulation_path)
    # get epsilon from HyperCLaw-V1.1 file [64, 64, length_z (but same values)]
    alpha = get_scalar_from_input_file('alpha', simulation_path)
    beta = get_scalar_from_input_file('beta', simulation_path)
    gamma = get_scalar_from_input_file('gamma', simulation_path)
    g11 = get_scalar_from_input_file('g11', simulation_path)
    material_properties = np.array([alpha, beta, gamma, g11])
    
    x_linarray = np.linspace(xyz_lo[0], xyz_hi[0], xyz_n_cell[0])
    y_linarray = np.linspace(xyz_lo[1], xyz_hi[1], xyz_n_cell[1])
    
    length_x = xyz_hi[0] - xyz_lo[0]
    length_y = xyz_hi[1] - xyz_lo[1]
    length_z = xyz_hi[2] - xyz_lo[2]
    nlayers = xyz_n_cell[2]
    thickness_per_layer = length_z / nlayers
    FE_lo, FE_hi = get_geometry_settings(simulation_path, component='FE', input_filename='inputc')
    FE_loz = FE_lo[2]
    FE_thickness = FE_hi[2] - FE_lo[2]
    index_fede_hi = 0
    index_fede_lo = 0
    for layer_idx in range(nlayers): # not super robust
        if np.round(layer_idx * thickness_per_layer, 11) > FE_loz:
            index_fede_hi = layer_idx
            index_fede_lo = layer_idx - 1
            #print(index_fede_hi*thickness_per_layer, index_fede_lo*thickness_per_layer, FE_loz)
            break
    #print(length_z, nlayers, index_fede_lo, index_fede_hi)

    SC_lo, SC_hi = get_geometry_settings(simulation_path, component='SC', input_filename='inputc')
    SC_hiz = SC_hi[2]
    #print(SC_hiz)
    index_sc_hi = 0
    #index_fede_lo = 0
    for layer_idx in range(nlayers): # not super robust
        if np.round(layer_idx * thickness_per_layer, 11) >= SC_hiz:
            index_sc_hi = layer_idx
            #index_fede_lo = layer_idx - 1
            #print(index_fede_hi*thickness_per_layer, index_fede_lo*thickness_per_layer, FE_loz)
            break
    
    #print(index_sc_hi)
    
    DE_lo, DE_hi = get_geometry_settings(simulation_path, component='DE', input_filename='inputc')
    DE_thickness = DE_hi[2] - DE_lo[2]
    
    #top_layer = xyz_n_cell[2] - 1

    voltage = initial_voltage # for python's sake
    bads = []
    goods = []
    for entry in os.scandir(simulation_path):
        if entry.is_dir() and 'plt' in entry.name:
            plt_name = entry.name
            ts = int(plt_name[3:])
            ds = yt.load(simulation_path + plt_name)
            ds.force_periodicity()
            ad = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
            try:
                P_array = ad['Pz'].to_ndarray()
                #polarization = P_array[:,:,top_layer]
                polarizations = P_array[:,:,:]
            except ValueError:
                print('ValueError in ', simulation_path, plt_name,'Pz. Skipping.')
                bads.append(plt_name)
                continue
            

            try:
                E_array = ad['Ez'].to_ndarray()
                #electric_field = E_array[:,:,top_layer]
                electric_fields = E_array[:,:,:]
            except ValueError:
                print('ValueError in ', simulation_path, plt_name,'Ez. Skipping.')
                bads.append(plt_name)
                continue
                
            try:
                Phi_array = ad['Phi'].to_ndarray()
                #electric_field = E_array[:,:,top_layer]
                #electric_fields = E_array[:,:,:]
            except ValueError:
                print('ValueError in ', simulation_path, plt_name,'Phi. Skipping.')
                bads.append(plt_name)
                continue
            
            try:
                charge_array = ad['charge'].to_ndarray()
                charge_map = charge_array[:,:,:]
            except ValueError:
                print('ValueError in ', simulation_path, plt_name,'Charges. Skipping.')
                bads.append(plt_name)
                continue
            
            #Calculate V_fe_avg
            V_FeDe = 0.5*(Phi_array[:,:,index_fede_lo] + Phi_array[:,:,index_fede_hi])
            # V_FeDe = Phi_array[:,:,index_fede_lo]
            integral_V = (1/length_x) * (1/length_y) * np.trapz(np.trapz(V_FeDe,x_linarray),y_linarray)
            
            #print(integral_V)
            
            
            Phi_sc_hiz = Phi_array[:,:,index_sc_hi]
            #print(np.min(Phi_sc_hiz.flatten()), np.max(Phi_sc_hiz.flatten()))
            #voltage_sc_hiz = (1/length_x) * (1/length_y) * np.trapz(np.trapz(Phi_sc_hiz,x_linarray),y_linarray)
            voltage_sc_hiz = np.array([np.min(Phi_sc_hiz.flatten()), np.max(Phi_sc_hiz.flatten())])
            
            epsilon_array = ad['epsilon'].to_ndarray()
            charges = np.zeros(nlayers) 
            for layer_idx in range(nlayers): # this whole thing only rigurously valid for top layer
                epsilon = epsilon_array[0,0,layer_idx] #should make this more robust
                E_field_contribution = epsilon*electric_fields[:,:,layer_idx]
                P_contribution = polarizations[:,:,layer_idx]
                D_field = E_field_contribution + P_contribution
                charges[layer_idx] = -1 * (1/length_x) * (1/length_y) *np.trapz(np.trapz(D_field,x_linarray),y_linarray)
            
            #print(len(voltage_changes_series.index))
            for r in range(1, len(voltage_changes_series.index)):
                if ts == voltage_changes_series.index[-1]:
                    voltage = voltage_changes_series.values[-1]
                    break
                if ts >= voltage_changes_series.index[r-1] and ts < voltage_changes_series.index[r]:
                    voltage = voltage_changes_series.values[r-1]
                    break
            #voltage = voltage_changes_series.values[-1]

          # if ts in voltage_changes_series.index:
          #     if ts == 0:
          #         pass
          #     else:
          #         voltage = voltage + voltage_step_increment
            
            filename = 'npz' + plt_name[3:] 
            
            thicknesses = np.array([10e-9, DE_thickness, FE_thickness])
            
            np.savez_compressed(simulation_path+filename, ts=ts, simulation_time=ts*dt, \
                        applied_voltage=voltage, integrated_charges=charges, \
                polarizations=polarizations, efields=electric_fields, Phi_SC=Phi_sc_hiz, \
                integrated_voltage=integral_V, switch_plts=voltage_changes_series.index, voltage_sc_hiz=voltage_sc_hiz, \
                               thicknesses=thicknesses, material_properties=material_properties, charge_map=charge_map)
            
            if verbose:
                print(simulation_path+filename)
                #print(ts, ts*dt, voltage, charge)#, polarization, electric_field)
                
            goods.append(plt_name)
    
    if len(bads) > 0:
        print('Processed', simulation_path, 'but with', len(bads), 'problems.')
    
    if len(goods) == 0:
        print('Nothing to process in', simulation_path, '.')
    else:
        print('Processed', len(goods),'plts in', simulation_path)
