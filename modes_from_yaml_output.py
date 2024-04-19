from pylab import *
import yaml
from yaml import Loader
import os

def getnatoms(poscar='POSCAR'):
    f = open(poscar,'r')
    lines = f.readlines()
    f.close()

    nat = int(sum(asfarray(lines[6].split())))

    return nat


def createShortenedYaml(nat):
    firstnlines = 12 + 3*nat + 6 + nat*(3+nat*4)
    #print('# Creating mesh-gamma.yaml')
    #print('# Skipping first %i lines'%firstnlines)

    os.system('head -n %i mesh.yaml > mesh-gamma.yaml'%firstnlines)

    return 

def getGammaEigenvectors(bands=array([3,4,5])):
    """
    prints eigenvectors for band indices in 'bands' array
    """

    stream = open("mesh-gamma.yaml", 'r')
    data = yaml.load(stream,Loader)
    natom = int(data['natom'])

    lattice = array(data['lattice'])
    points = array(data['points'])

    eigenvectors = zeros([size(bands),natom,3])
    
    for j, band in enumerate(bands):
        #print('# Eigenvector for band index = %i'%band)
        for i in range(natom):
            mode = asarray(data['phonon'][0]['band'][band]['eigenvector'][i])[:,0]
            eigenvectors[j,i,:] = mode
            #print("%16.12f %16.12f %16.12f"%(mode[0],mode[1],mode[2]))
            
    return eigenvectors

def writeVestaMode(bands, eigvec, nat, scaling_factor):

    sf = scaling_factor


    for j, band in enumerate(bands):
        
        towrite = "VECTR\n"
        for i in range(1,1+nat):
            towrite += "%4d%9.5f%9.5f%9.5f 0\n"%(i,eigvec[j][i-1][0]*sf,eigvec[j][i-1][1]*sf,eigvec[j][i-1][2]*sf)
            towrite += "%5d  0   0    0    0\n 0 0 0 0 0\n"%i
        towrite += " 0 0 0 0 0\n"
        towrite += "VECTT\n"
        for i in range(1,1+nat):
            towrite += "%4d%6.3f 255   0   0 1\n"%(i,0.5)
            towrite += " 0 0 0 0 0\n"


        print('------------------------')
        print('mode %i'%band)
        print('------------------------')
        print(towrite)
    

if __name__ == '__main__':

    nat = getnatoms('POSCAR')
    createShortenedYaml(nat)

    bands = array([3,4,5])
    evecs = getGammaEigenvectors(bands)

    sf = 4
    writeVestaMode(bands,evecs,nat,sf)
