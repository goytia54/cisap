from __future__ import division
import hoomd
from hoomd import md
from hoomd import deprecated
import gsd.hoomd
import settings
import random
import numpy as np

'''####################################################
# hoomdRun() creates a HOOMD input file and then calls 
# HOOMD to run MD simulations. Parameters are explained
# below
####################################################'''


def hoomdRun(A,B,C,Idiag,tM,run_type):
	tempi = settings.temp #Inital temp of run
	thermo = settings.thermo #how often to output the energies
	movie = settings.dump #how often to output for visualization
	steps = settings.steps # number of runs at each temp
	dT = settings.tstep #time step of simulation
	eps = settings.eps #epsilons
	sig = settings.sig #sigmas
	rS = settings.rS # how often to write restart
	rc = max(sig)*2.5 # max cut off for sigmas (2.5*largest sigma)
	if run_type == -1: #mix dimer initialization
		rSname = 'run_mix.rS'
		log_name = 'random_mix.log'
		initfile = 'init_mix.gsd'
		finalxml = 'random_mix.xml'
		finaldcd = 'random_mix.dcd'
		flip = -1 # flips the coordinates to make other handedness
		D_type = '5' #add another type of middle atom
		P_id = 1 # particle ID
		ijrange = 5 #number of differnt atoms
		types = ['1','2','3','4','5']
	else:   #pure dimer intialization
		rSname = 'run_pure.rS'
		log_name = 'random_pure.log'
		initfile = 'init_pure.gsd'
		finalxml = 'random_pure.xml'
		finaldcd = 'random_pure.dcd'
		flip = 1
		D_type = '4'
		P_id = 0
		ijrange = 4
		types = ['1','2','3','4']
'''###############################################################
#  Creation of the input file for the HOOMD run using .gsd format
#  Further documentation can be found on HOOMD page
##############################################################'''	

	s = gsd.hoomd.Snapshot() 
	s.particles.N = 2 # center of mass particles
	numP = s.particles.N #number of particles
	s.particles.types = ['4',D_type] #intialize particles in which dtype can be '4' or '5'
	s.particles.typeid = [0,P_id]
	s.particles.position = [[-3.5,0,0],[3.5,0,0]]
	s.particles.body = [0,1]
	s.particles.mass = 2*[tM]
	s.particles.charge = [0,0]
	s.particles.diamter = [0.05,0.05]
	s.particles.moment_inertia = [numP*[Idiag]]
	s.configuration.box = [14,14,14,0,0,0]
	s.configuration.dimensions = 2
	gsd.hoomd.create(name=initfile,snapshot = s)
	

'''#############################################################
#  Begin running the MD simulation
#
############################################################'''

	hoomd.context.initialize("")
	system = hoomd.init.read_gsd(filename=initfile) #reads in the file
	rcut =  max(sig)*2.5
	nl = md.nlist.cell(check_period=1) # neighbor list
	rigid = hoomd.md.constrain.rigid() #creates rigid group
	# adds the surround atoms  
	system.particles.types.add('1') 				
	system.particles.types.add('2') 
	system.particles.types.add('3') 
	# assign the atoms to rigid COM atom 
	rigid.set_param('4',positions=[(A[0],A[1],A[2]),(B[0],B[1],B[2]),(C[0],C[1],C[2])],types=['1','2','3'],diameters=[sig[0],sig[1],sig[2]])
	# adds the fith type if using the mixed dimmer and flips the coordinates
	if run_type == -1:
		rigid.set_param('5',positions=[(flip*A[0],A[1],A[2]),(flip*B[0],B[1],B[2]),(flip*C[0],C[1],C[2])],types=['1','2','3'],diameters=[sig[0],sig[1],sig[2]])
	rigid.enable() 
	rigid.create_bodies() #creates rigid body
	center = hoomd.group.rigid_center()# used for rigid body integration
	part = hoomd.group.all() #used for dumping all particles
	lj = md.pair.lj(r_cut=rcut,nlist=nl) #Leonard-Jones parameters set up
	count = 0
	# setting all the LJ interactions for number of particles 
	for i in range(0,ijrange):
		for j in range(i,ijrange):
				if i == 3 or j == 3 or i == 4 or j == 4:
					# interactions wih type 4 and 5 or COM is 0
					lj.pair_coeff.set(types[i],types[j],epsilon=0,sigma=0.1) 
				else:
					lj.pair_coeff.set(types[i],types[j],epsilon=0.5*(eps[i]+eps[j]),sigma=0.5*(sig[i]+sig[j]))
					count += 1	
	hoomd.md.update.enforce2d() #enforcing 2D dynamics
	hoomd.md.integrate.mode_standard(dt=dT) # setting timestep for integrator
	rNUM = random.randint(1,100000) #random number generator for seed
	lang = hoomd.md.integrate.langevin(group=center, kT=tempi, seed = rNUM) # langevin dynamics
        
	# outputs, check HOOMD documentation for further explanation
	hoomd.analyze.log(filename=log_name,quantities=['temperature','kinetic_energy','potential_energy','pressure','volume'],period=thermo,header_prefix='#',overwrite=True)
	hoomd.deprecated.dump.xml(group=part,filename =finalxml,vis=True,image=True)#period=100)
	hoomd.dump.dcd(filename=finaldcd,period=movie, group = part ,overwrite = True)
	
	# This is the method to change the temperature every so many steps as well as 
	# decreasing the rate at which the temperature changes to cool slower and slower.
	hoomd.run(100000)
	rate = 0.90
	while steps > 0: #number of steps at each TEMP, this is in the settings.py file
		tempi = rate*tempi
		lang.set_params(kT=tempi)
		hoomd.run(100000)
		steps -= 1
		rate += 0.0012
	lang.set_params(kT=0.001)
	hoomd.run(100001)


