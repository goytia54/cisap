import settings
import random

'''#############################################
#   EorS is a function used to change the 
#   parameters of the trimer molecules to determine
#   what parameters tend to cause the lower
#   energy of like dimers.
###########################################'''


def EorS(): # returns sigmas

#  Intialize  parameters 
    eps = settings.eps
    sig = settings.sig
    length = settings.length
    r_angle = settings.angle
    dSig = settings.dSig
    dEps = settings.dEps
    dA = settings.dAng
    dL = settings.dLen
    
'''########################################################################
Determines which parameter to change, varry the randint range if you include
further parameters
#########################################################################'''
    rand_parm = random.randint(0,1) 
    
    skip = 0

'''########################################################################## 
Rand_p will need to be adjusted from the randint as it describes the index
of what you wil be changing. 
	-i.e. is randit is from (0,4) and 2,3,4 correspond to changing the sigmas
	then randp will need to have 2 taken away from randit in order to have the
	proper index
##############################################################################'''
    
    rand_p = rand_parm
    
    delt_eps = dEps*(2*random.random()-1)  # variable to either increase or decrease paramter
    mol_eps = delt_eps + eps[rand_p]  # add the change
    if mol_eps > 0 and mol_eps <= 1.0: # determine if change is between range set
    	settings.eps[rand_p] = mol_eps # if it works then change is accepted and HOOMD is run
        skip = 0 
    else: #else HOOMD is not ran and last energy is recorded
	skip = 1

############################################
# Uncomment if you want to change the sigmas, 
# length between atoms, and/or angle between
# atoms.
#


#                 SIGMAS
#    elif rand_parm >= 0 and rand_parm < 3:
#        rand_p = rand_parm
#        delt_sig = dSig*(2*random.random()-1)
#        mol_sig = delt_sig+ sig[rand_p]
#        if mol_sig >= 0.5 and mol_sig <= 2.0:
#        	settings.sig[rand_p] = mol_sig
#		skip = 0
#        else:
#		skip = 1

#           	 LENGTHS
#    elif rand_parm == 9 or rand_parm == 10: #change length between
#        rand_p = rand_parm - 9
#        delt_len = dL*(2*random.random()-1)
#        new_l = delt_len + length[rand_p]
#        if new_l >= 0.5 and new_l <= 3.0:
#                settings.length[rand_p] = new_l
#                skip = 0
#        else:
#                skip = 1

#		ANGLES
#    elif rand_parm == 11: #change angle
#        delt_ang = dA*(2*random.random()-1)
#        r_angle  = r_angle + delt_ang
#        if r_angle >= 0 and r_angle <= 180:
#                settings.angle =  r_angle
#                skip = 0
#        else:
#                skip = 1
    return skip #returns if move allowed or not

