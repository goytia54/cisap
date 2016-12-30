import settings
import math as m
import random
import numpy as np

'''########################################################
# mol_xyz() returns the coordiantes of the molecule and 
# accounts for changes in all paramters(sigmas,length, and 
# angle)
########################################################'''

def mol_xyz():
    #intialize parametrs
    sig = settings.sig
    siga = sig[0]/2
    sigb = sig[1]/2
    sigc = sig[2]/2
    dA = settings.dAng
    dL = settings.dLen
    length = settings.length
    r_angle = settings.angle
    
    # adjust angles
    if r_angle <= 90:
        rand_rad = r_angle/57.2958
	flip = -1
    else:
        rand_rad = (180-r_angle)/57.2958
	flip = 1
    
    # a,b,c are the lengths in the triangle where a=CB,b=AC,c=AB
    # we don't know the distance between AC or b because angle
    # and lengths can change so need to calculate and then find the
    # coordinates of C bases on A and B
    
    c = settings.length[0]
    a = settings.length[1]
    b = (a**2)+(c**2)-2*a*c*m.cos(rand_rad)
    Cy = a*m.sin(rand_rad)
    Cx = flip*m.sqrt(a**2-Cy**2)
    A = [-c,0,0]
    B = [0,0,0]
    C = [Cx,Cy,0]
    return A, B, C

'''##########################################################
#  com() taks in the postions of the three Atoms and determines
#  the principle momemnt of inertia from diagonalizing 
#  inertia matrix in page 164 of "Spectra of Atoms and Molecules" 
#  2nd Ed. by Peter F. Bernath
##########################################################'''


def com(A,B,C,flip):
	M = []
	TotalM = 0
	sig = settings.sig
	den = 1.40256
	for i in range(0,3):
	        mass = den*(4/3)*(sig[i]**3)*3.141592653589793 #determine all the masses of the atoms
		M.append(mass) 
		TotalM += mass

	#determining the center of mass for the molecules
	cX = (M[0]*A[0] + M[1]*B[0] + M[2]*C[0])/TotalM #COM coord, X
	cY = (M[0]*A[1] + M[1]*B[1] + M[2]*C[1])/TotalM # Y
	cZ = (M[0]*A[2] + M[1]*B[2] + M[2]*C[2])/TotalM # Z
	com = [cX,cY,cZ]
	cXYZ = [[A[0]-cX,A[1]-cY,A[2]-cZ],[B[0]-cX,B[1]-cY,B[2]-cZ],[C[0]-cX,C[1]-cY,C[2]-cZ]] # new coords
	Ixx,Ixy,Ixz,Iyy,Iyz,Izz = 0.0,0.0,0.0,0.0,0.0,0.0

	#principle points of Inertia matrix, see book for further detials
	for i in xrange(3):
		Ixx += M[i]*(cXYZ[i][1]**2 + cXYZ[i][2]**2)
		Iyy += M[i]*(cXYZ[i][0]**2 + cXYZ[i][2]**2)
		Izz += M[i]*(cXYZ[i][0]**2 + cXYZ[i][1]**2)
		Ixy += -M[i]*cXYZ[i][0]*cXYZ[i][1]
		Ixz += -M[i]*cXYZ[i][0]*cXYZ[i][2]
		Iyz += -M[i]*cXYZ[i][1]*cXYZ[i][2]

	#construct matrix
	Imatrix = np.matrix([[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]])
	
	#diagonalize the matrix
	Idiag,imatrix = np.linalg.eig(Imatrix)
	comX = 0-com[0]
	comY = 0-com[1]

	#returns new coordiantes of atoms in relation to center of mass, needed for HOOMD
	# initialization
	A = [A[0]+comX,A[1]+comY,A[2]]
	B = [B[0]+comX,B[1]+comY,B[2]]
	C = [C[0]+comX,C[1]+comY,C[2]]
	return A,B,C,Idiag,TotalM
