#!/usr/bin/env python3
#! -*- coding: utf-8 -*-

"""
Plot TTT (or CCT) diagram and transformed fraction
"""
import matplotlib.pyplot as plt
from transformation_models import Alloy, TransformationDiagrams
import numpy as np
from math import log, exp, expm1
import pandas as pd


if __name__ == '__main__':
    # Defines alloy (grain size gs and composition)
    alloy = Alloy(gs=3, C=0.35, Mn=0.37)
    # alloy = Alloy(gs=7, C=1, Mn=0.7, Si=0.3, Ni=0.15, Cr=5.125, Mo=1.15)

    # Initializes diagrams object
    diagrams = TransformationDiagrams(alloy)

    """    
    Example 1: Plot TTT and CCT diagrams
    """
    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig1.subplots_adjust(wspace=.2)

    diagrams.TTT(ax=ax1)
    #ax1.set_xlim(1e-2, 1e8)
    #ax1.set_ylim(300, 1000)

    diagrams.CCT(ax=ax2)
    # diagrams.CCT(ax=ax2, phi_min=900/1e19, phi_max=900/1e8)
    #ax2.set_xlim(1e-2, 1e8)
    #ax2.set_ylim(300, 1000)

    fig1.suptitle(ax1.get_title())
    ax1.set_title('')
    ax2.set_title('')

    """
    Example 2: Get phase fraction information and save data in file
    """
    # Chooses a thermal cycle (continuous cooling from 1000 to 0 oC in
    # 2000 s), draws cooling curve in the CCT and plots phase fraction
    t, T = [0, 10e4], [843, 0]

    # get_transformed_fraction returns a pandas DataFrame
    # n is number of points
    data = diagrams.get_transformed_fraction(t, T, n=2001)

    # Displays data
    print(data)

    # Save as csv
    data.to_csv('data.csv', index=False)

    # Save as excel
    #data.to_excel('data.xls', index=False)
    tempiObbiettivo = diagrams.df_ferrite

    ts_ferrite = tempiObbiettivo.iloc[:, 1]
    tf_ferrite = tempiObbiettivo.iloc[:, 2]

    tau_ferrite = (log(0.01/0.9)-tf_ferrite-ts_ferrite)/log(1-0.01/0.9)
    xsi_eq_ferrite = 0.01/(1-np.exp(-np.divide(ts_ferrite, tau_ferrite)))
    # xsi_eq_ferrite = np.ones(len(tau_ferrite))

    tm_ferrite = -tau_ferrite*np.log(1-np.divide(0.5, xsi_eq_ferrite))

    tempiObbiettivo.insert(3, "tm_ferrite", tm_ferrite, True)
    tempiObbiettivo.insert(4, "tau_ferrite", tau_ferrite, True)
    tempiObbiettivo.insert(5, "xsi_eq_ferrite", xsi_eq_ferrite, True)

    print(tempiObbiettivo)

    tempiObbiettivo.to_csv('ferrite.csv', index=False)

    tempiObbiettivo = diagrams.df_perlite

    ts_perlite = tempiObbiettivo.iloc[:, 1]
    tf_perlite = tempiObbiettivo.iloc[:, 2]

    tau_perlite = (log(0.01/0.9)-tf_perlite-ts_perlite)/log(1-0.01/0.9)
    xsi_eq_perlite = 0.01/(1-np.exp(-np.divide(ts_perlite, tau_perlite)))
    # xsi_eq_perlite = np.ones(len(tau_perlite))

    tm_perlite = -tau_perlite*np.log(1-np.divide(0.5, xsi_eq_perlite))

    tempiObbiettivo.insert(3, "tm_perlite", tm_perlite, True)
    tempiObbiettivo.insert(4, "tau_perlite", tau_perlite, True)
    tempiObbiettivo.insert(5, "xsi_eq_perlite", xsi_eq_perlite, True)

    print(tempiObbiettivo)

    tempiObbiettivo.to_csv('perlite.csv', index=False)

    tempiObbiettivo = diagrams.df_bainite

    ts_bainite = tempiObbiettivo.iloc[:, 1]
    tf_bainite = tempiObbiettivo.iloc[:, 2]

    tau_bainite = (log(0.01/0.9)-tf_bainite-ts_bainite)/log(1-0.01/0.9)
    xsi_eq_bainite = 0.01/(1-np.exp(-np.divide(ts_bainite, tau_bainite)))
    # xsi_eq_bainite = np.ones(len(tau_bainite))

    tm_bainite = -tau_bainite*np.log(1-np.divide(0.5, xsi_eq_bainite))

    tempiObbiettivo.insert(3, "tm_bainite", tm_bainite, True)
    tempiObbiettivo.insert(4, "tau_bainite", tau_bainite, True)
    tempiObbiettivo.insert(5, "xsi_eq_bainite", xsi_eq_bainite, True)

    print(tempiObbiettivo)

    tempiObbiettivo.to_csv('bainite.csv', index=False)

    """
    Example 3: Plot CCT diagram and transformed fraction
    """
    fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 6))
    fig2.subplots_adjust(wspace=.2)

    # Plot CCT diagram
    diagrams.CCT(Tini=1000, ax=ax3)

    # Same thermal cycle as from example 2
    diagrams.draw_thermal_cycle(ax3, t, T)
    diagrams.plot_phase_fraction(t, T, xaxis='T', ax=ax4)

    fig2.suptitle(ax3.get_title())
    ax3.set_title('')
    ax4.set_title('')

    plt.show()
