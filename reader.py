import settings

'''####################################################
#   reader.py
#  -Reads to the HOOMD output files and then extracts 
#   the potenatial energy and returns it to be read 
#   into SimAneal for comparison
####################################################'''    

def readTraj(filename):
	with open(filename,'rb') as f:
    		f.seek(-2, 2)             # Jump to the second last byte.
    		while f.read(1) != b"\n": # Until EOL is found...
        		f.seek(-2, 1)         # ...jump back the read byte plus one more.
    		last = f.readline() # splits up data
		last = last.split()
        f.close()
	return last[3]

