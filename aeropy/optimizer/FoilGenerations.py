import os
from aeropy import xfoil_module as xf
import numpy as np
from shutil import copy2

class flapFoil:
    projectName="projectNameNotSet"
    generationName ="generationNameNotSet"
    rawFoilName = "rawFoilNameNotSet"
    hingedFoilName = "hingedFoilNameNotSet"
    flapFoilName = "flapFoilNameNotSet"
    
    #Sets the relative path of the flapped Foil
    """ def getRelFoilPath(self):
        relFoilPath = "airfoilGens" + '/' + generationName + '/' + rawFoilName + '/' + hingedFoilName + '/' + flapFoilName + '.dat'
        return relFoilPath """


    #Constructor for class flapFoil
    def __init__(self, inputProjectName, inputGenerationName, inputRawFoilName, inputHingedFoilName,  inputFlapFoilName):
        self.flapFoilName = inputFlapFoilName
        self.hingedFoilName = inputHingedFoilName
        self.rawFoilName = inputRawFoilName
        self.generationName = inputGenerationName
        self.projectName = inputProjectName


        if os.path.isdir(self.flapFoilDir()):
            pass
            #print("The directory of the flapped Foil is ok!")

        else:
            pass
            #print("The directory of the flapped foil does not exist")

        
            try:
                os.makedirs(self.flapFoilDir())
            except OSError:
                print ("Error creating directory %s failed" % self.flapFoilDir())
        
        #print("Created flapFoil %s.%s.%s.%s.%s" %(self.projectName, self.generationName, self.rawFoilName, self.hingedFoilName, self.flapFoilName))

    def flapFoilDir(self):
        return os.getcwd() + "/airfoilGens/" + self.generationName + '/' + self.rawFoilName +  '/' + self.hingedFoilName
        
    def flapFoilPath(self):
        return self.flapFoilDir() +  '/' + self.flapFoilName

    def relFlapFoilPath(self):
        return "airfoilGens/" + self.generationName + '/' + self.rawFoilName + '/' + self.hingedFoilName + '/' + self.flapFoilName

    def relFlapFoilDir(self):
        return "airfoilGens/" + self.generationName + '/' + self.rawFoilName + '/' + self.hingedFoilName + '/'

    def analyzeFoil(self, inputAlpha, inputReynolds, inputIteration):
        if type(inputReynolds) == list or type(inputReynolds) == np.ndarray:
            for re in inputReynolds:
                xf.find_coefficients(airfoil=self.flapFoilName, alpha=inputAlpha, Reynolds=re, iteration=inputIteration, \
                    NACA=False, delete=False, dir=self.relFlapFoilDir() + '/')
        else:
            xf.find_coefficients(airfoil=self.flapFoilName, alpha=inputAlpha, Reynolds=inputReynolds, iteration=inputIteration, \
                NACA=False, delete=False, dir=self.relFlapFoilDir())

    def loadFoil(self, inputFoilToLoad):
        copy2(os.getcwd() + '/foilsToLoad/' + inputFoilToLoad, self.flapFoilPath() )




class hingedFoil:
    projectName="projectNameNotSet"
    generationName ="generationNameNotSet"
    rawFoilName = "rawFoilNameNotSet"
    hingedFoilName = "hingedFoilNameNotSet"

class rawFoil:
    projectName="projectNameNotSet"
    generationName ="generationNameNotSet"
    rawFoilName = "rawFoilNameNotSet"

class generation:
    projectName="projectNameNotSet"
    generationName ="generationNameNotSet"


def main():
    airf1 = flapFoil('EPA_0', 'G0', 'R0', 'H0', 'F0' )
    #print(airf1.flapFoilPath())
    airf1.loadFoil("FE_H5.dat")
    reArray=np.arange(1e5, 3e6, 1e5)
    alphaArray=np.arange(0, 30, 1)
    airf1.analyzeFoil(alphaArray, reArray, 20)
    
if __name__ == "__main__":
    main()
