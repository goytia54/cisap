import math as m
import random
import settings
import os

def simAneal(eng_mix,eng_pure,skip):
	oldDE = settings.old_dEng
	newDE = float(eng_pure)-float(eng_mix)
        kT_var = settings.kT
	Expon = m.exp(-(newDE-oldDE)/kT_var)
	eps = settings.eps
	oldEps = settings.old_eps
	sig = settings.sig
	oldSig = settings.old_sig
	ang = settings.angle
	lenth_n = settings.length
	oldang = settings.old_angle
	oldlen = settings.old_length
	rand_parm = random.random()
	if newDE < oldDE:
                settings.old_sig = [i for i in settings.sig]
                settings.old_eps = [i for i in settings.eps]
                settings.old_dEng = newDE
                settings.old_length = [i for i in settings.length]
                settings.old_angle = ang
                settings.accept += 1
        elif Expon > rand_parm:
                settings.old_sig = [i for i in settings.sig]
                settings.old_eps = [i for i in settings.eps]
                settings.old_dEng = newDE
                settings.old_length = [i for i in settings.length]
                settings.old_angle = ang
                settings.accept += 1
        else:
                settings.sig = [i for i in settings.old_sig]
                settings.eps = [i for i in settings.old_eps]
                settings.length = [i for i in settings.old_length]
                settings.angle = oldang
	if settings.rounds == 0:
        	fit_file = open('fit.txt','w')
        	fit_file.write('# rounds eps1 eps2 ep3 newDE oldDE kT_var accept skip \n')
		fit_file.write('#    1     2   3    4   5     6      7      8      9   \n') 
    	else:
        	fit_file = open('fit.txt','a')
	oldEps = settings.old_eps
        oldSig = settings.old_sig
    	fit_file.write('{0} '.format(settings.rounds))
	fit_file.write('{0} {1} {2} '.format(oldEps[0],oldEps[1],oldEps[2]))
	fit_file.write('{0} {1} {2} {3} {4}\n'.format(newDE,oldDE,kT_var,settings.accept,skip))
    	fit_file.close()
	if not os.path.exists("good_traj/"):
        	os.system('mkdir good_traj')
    	os.chdir('good_traj')
	currentRound = settings.rounds
        os.system('cp ../random_mix.xml random_mix{0}.xml'.format(currentRound))
        os.system('cp ../random_mix.dcd random_mix{0}.dcd'.format(currentRound))
 	os.system('cp ../random_pure.xml random_pure{0}.xml'.format(currentRound))
        os.system('cp ../random_pure.dcd random_pure{0}.dcd'.format(currentRound))
	os.system('cp ../random_mix.log random_mix{0}.log'.format(currentRound))
        os.system('cp ../random_pure.log random_pure{0}.log'.format(currentRound))

	os.chdir('..')
