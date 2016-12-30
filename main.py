# -*- coding: utf-8 -*-
"""
Created on Mon May 23 16:17:44 2016

@author: Michael Goytia 
"""

'''##############################################################
#            CHIRAL DIMER SIMULATED ANNEALING PROCEDURE (CISAP) 
#     		-A Monte Carlo-type alogrithim
# 1. Varry the paramters of a trimer (sigmas, epsilons) randomly
#    and determin if the variation is accepted or not
#    A. If the move is accepted then both handedness input files
# 	are constructed to be run in HOOMD (skip == 0) 
#    B. If the move is rejected then record the energy of the 
#	the last trajectory
# 2. Run the HOOMD simulations twice for each configuration,
#    Pure and Mixed dimers and then determine if energies 
#    differences are better than the previous run.
#    A. If energies are BETTER than the last run of fall between
#	the allowed difference based on evolutionary temperature
#	then new parameters are accpeted
#    B. If energies are WORST than the previous run then 
#       then parameters are not accepted and return to step one
# 3. Energies(new old), paramters, skipped runs, and accepted run 
#    are found in fit.txt
# 4. Evolutionary temperature(kT) decreases over epochs (~200-300 epochs)
#
#  !!! RUNS currently wiht HOOMD 2.X.X and to have HOOMD module loaded 
#  from the server in order to run
#
#  !!! run "python main.py" 
####################################################################'''   

from __future__ import division
import random 
import math as m
import os
from launch import *
from molecalc import *
from reader import *
from settings import *
from sim_aneal import *
from run2 import *
import subprocess

def main():
	rounds = settings.rounds #rounds
    	max_inter = settings.max_inter
   	while rounds < max_inter:
		kT = settings.kT #Evolutionary temperature
        	if rounds % 200 == 0 and rounds > 0: #decrease Evolutionary temperature every 200 steps
        		settings.kT =  kT * 0.9
		skip = EorS() #pick a random paramter to change then returns if move is accepted
		if skip == 0: #runs if parameter chnage is accpeted 
			A,B,C = mol_xyz() # returns A,B,C,D
			A1,B1,C1,Idiag,TotalM = com(A,B,C,1) #determines the moment of inertia 
			en_pure = []
        	        en_mix = []
			for i in range(2):
				hoomdRun(A1,B1,C1,Idiag,TotalM,-1) #run hoomd for mix
				hoomdRun(A1,B1,C1,Idiag,TotalM,1)  #run hoomd for the pure
				energy_mix = readTraj('random_mix.log') #extract energires for comparison
				energy_pure = readTraj('random_pure.log')
				en_pure.append(float(energy_pure)) #gathering two runs for each dimer
                	        en_mix.append(float(energy_mix))
			energy_pure,energy_mix=min(en_pure),min(en_mix) #takes the best or lowest energy dimer
			simAneal(energy_mix,energy_pure,skip) # determine if the difference in energy is better
		else: #if the move is not accepted from EorS()
			energy_pure,energy_mix=min(en_pure),min(en_mix) 
			simAneal(energy_mix,energy_pure,skip)	
        	settings.rounds += 1 #keep tracks of the number of EPOCHS
        	rounds = settings.rounds

if __name__=="__main__":
     main()

