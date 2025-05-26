########################################################################################################
#                                                                                                      #
#  Project:  TFG chemistry                                                                             #
#  Author:   Domènec Huerta Estradé                                                                    #
#  Date:     03/01/2025                                                                                #
#  Purpose:  Definig useful functions used generally in all the program like plotting graphs           #
#                                                                                                      #
########################################################################################################

import os
import numpy as np
import matplotlib.pyplot as plt

def folder_check(folder_name):
    """
    Checks for a required folder required in the code ans
    """
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

def ploter(freq_1, freq_2, directory, name, save_format, ploter_show = False):
    if len(freq_2) != len(freq_1):
        for i in range(len(freq_2) - len(freq_1)):
            freq_1.append(0)
    
    diference_vec = [freq_2[i] - freq_1[i] for i in range(len(freq_1))]
    print(freq_1)
    print(diference_vec)
    # Bins for both histograms
    bins = np.arange(1,len (freq_1)+1)
    
    # Create figure
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    # Histograms
    ax1.bar(bins, freq_1, alpha=1, width =1, color='blue', edgecolor='black', label="Database frequencies")
    ax1.set_xlabel("Conformation")
    ax1.set_ylabel("Database frequencies", color='blue')
    ax1.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax1.set_yscale("log")
    # Create second axis sharing the same x-axis
    ax2 = ax1.twinx()
    
    # Second histogram (Logarithmic scale)
    ax2.bar(bins,diference_vec, alpha=0.6, color='red', edgecolor='black', label="Analized protein frequencies")
    ax2.set_ylabel("Analized protein frequencies", color='red')
    ax2.set_yscale("log")
    # Legends for both datasets
    fig.legend(loc="upper right", bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)

    if not os.path.exists(directory):
        os.makedirs(directory)
    plt.savefig(directory+"/"+name+"."+save_format, format=save_format)
    if (ploter_show ==True):
        plt.show()
    max_len = len(freq_1)
    plt.xticks(range(1, max_len, (max_len//50)+1))
    plt.close()

def triple_ploter(freq_1, freq_2, freq_3, directory, name, save_format, ploter_show = False, masking = False, ordering = False, title = True):
    #Same function as plotter 
    max_len = max(len(freq_3), len(freq_2))
    
    if max_len != len(freq_1):
        for i in range(max_len - len(freq_1)):
            freq_1.append(0)
    
    if max_len != len(freq_2):
        for i in range(max_len - len(freq_2)):
            freq_2.append(0)
            
    if max_len != len(freq_3):
        for i in range(max_len - len(freq_3)):
            freq_3.append(0)
    
    diference_1 = [freq_2[i] - freq_1[i] for i in range(len(freq_1))]
    diference_2 = [freq_3[i] - freq_1[i] for i in range(len(freq_1))]
    
    if masking == True:
        mask = np.logical_not((np.array(diference_1) == 0) & (np.array(diference_2) == 0))
        freq_1 = np.array(freq_1)[mask]
        diference_1 = np.array(diference_1)[mask]
        diference_2 = np.array(diference_2)[mask]
    print(freq_1, diference_1, diference_2)
    if ordering == True:
        sorted_indices = np.argsort(freq_1)[::-1]
        freq_1 = np.array(freq_1)[sorted_indices]
        diference_1 = np.array(diference_1)[sorted_indices]
        diference_2 = np.array(diference_2)[sorted_indices]
    # Bins for both histograms
    bins = np.arange(1,len(freq_1)+1)
    bins_1 = [k - 0.2 for k in bins]
    bins_2 = [k + 0.2 for k in bins]    
    # Create figure
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.grid(True, which="both", linestyle="--", linewidth=0.4)
    # Histograms
    ax1.bar(bins, freq_1, alpha=0.9, width =0.9, color='#6c8ebf', edgecolor='black', label="Database conformation frequency")
    ax1.set_xlabel("Conformation")
    ax1.set_ylabel("Database conformation frequency", color='black')

    ax1.set_yscale("log")
    # Create second axis sharing the same x-axis
    ax2 = ax1.twinx()
    plt.xticks(range(1, max_len, int (max_len/40 )+1))
    # Second histogram (Logarithmic scale)
    ax2.bar(bins_1,diference_1, alpha=0.9, width =0.4, color='#daa520', edgecolor='black', label="Experimental conformation frequency")
    ax2.bar(bins_2,diference_2, alpha=0.9, width =0.4, color='#3cb371', edgecolor='black', label="Modeled conformation frequency")
    ax2.set_ylabel("Analyzed structure conformation frequency", color='black')
    ax2.set_yscale("log")
    # Legends for both datasets
    fig.legend(loc="upper right", bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)
    if (title == True):
        plt.title("Freqüencies de conformacions per a l'aminoàcid "+name+ " d'entre les "+str(max_len)+" a la base de dades")
    if (title == False):
        print("Freqüencies de conformacions per a l'aminoàcid "+name+ " d'entre les "+str(max_len)+" a la base de dades")
    if not os.path.exists(directory):
        os.makedirs(directory)
    plt.savefig(directory+"/"+name+"."+save_format, format=save_format)
    if (ploter_show ==True):
        plt.show()
    plt.close()

