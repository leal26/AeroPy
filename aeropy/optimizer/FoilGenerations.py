import os
from aeropy import xfoil_module as xf
import numpy as np
from shutil import copy2


class generation:
    projectName = "projectNameNotSet"
    generationName = "generationNameNotSet"

    def __init__(
        self,
        projectName,
        generationName):
        self.projectName = projectName
        self.generationName = generationName

        if os.path.isdir(self.generationDir()):
            pass
            print("The directory of the foil generation is ok!")

        else:
            # print("The directory of the foil generation does not exist")
            try:
                os.makedirs(self.generationDir())
            except OSError:
                print(
                    "Error "
                )
        


    def generationDir(self):
        return os.getcwd() + "/airfoilGens/" + self.generationName + "/"

    def relGenerationDir(self):
        return "airfoilGens/" + self.generationName + '/'





class rawFoil(generation):
    rawFoilName="rawFoilNameNotSet"

    def __init__(
        self,
        rawFoilName):
        self.rawFoilName=rawFoilName

        if os.path.isdir(self.rawFoilDir()):
            pass
            # print("The directory of the raw Foil is ok!")

        else:
            # print("The directory of the raw foil does not exist")
            try:
                os.makedirs(self.rawFoilDir())
            except OSError:
                print(
                    "Error "
                )
        


    def rawFoilDir(self):
        return os.getcwd() + "/airfoilGens/" + self.generationName + "/" + self.rawFoilName + "/"

    def relRawFoilDir(self):
        return "airfoilGens/" + self.generationName + '/' + self.rawFoilName + "/"

    def loadFoil(self, foilToLoad):
        copy2(
            os.getcwd() +
            '/foilsToLoad/' +
            foilToLoad,
            self.rawFoilDir())

class hingedFoil(rawFoil):
    hingedFoilName = "hingedFoilNameNotSet"

    def __init__(
        self,
        hingedFoilName):
        self.hingedFoilName = hingedFoilName

        if os.path.isdir(self.hingedFoilDir()):
            pass
            # print("The directory of the hinged Foil is ok!")

        else:
            # print("The directory of the hinged foil does not exist")
            try:
                os.makedirs(self.hingedFoilDir())
            except OSError:
                print(
                    "Error "
                )

        def hingedFoilDir(self):
            return os.getcwd() + "/airfoilGens/" + self.generationName + "/" + self.rawFoilName + "/" + self.hingedFoilName + "/"

        def relHingedFoilDir(self):
            return "airfoilGens/" + self.generationName + '/' + self.rawFoilName + "/" + self.hingedFoilName + "/"
        

class flapFoil(hingedFoil):
    flapFoilName="flapFoilNameNotSet"


    # Sets the relative path of the flapped Foil
    """ def getRelFoilPath(self):
        relFoilPath = "airfoilGens" + '/' + generationName + '/'
        + rawFoilName + '/' + hingedFoilName + '/' + flapFoilName + '.dat'
        return relFoilPath """

    # Constructor for class flapFoil
    def __init__(
            self,
            flapFoilName):
        self.flapFoilName=flapFoilName

        if os.path.isdir(self.flapFoilDir()):
            pass
            # print("The directory of the flapped Foil is ok!")

        else:
            pass
            # print("The directory of the flapped foil does not exist")

            try:
                os.makedirs(self.flapFoilDir())
            except OSError:
                print(
                    "Error creating directory %s failed" %
                    self.flapFoilDir())


    def flapFoilDir(self):
        return self.flapFoilDir() + '/' + self.flapFoilName

    def relFlapFoilDir(self):
        return "airfoilGens/" + self.generationName + '/' + self.rawFoilName + \
            '/' + self.hingedFoilName + '/' + self.flapFoilName


    def analyzeFlapFoil(self, inputAlpha, inputReynolds, inputIteration):
        if isinstance(
                inputReynolds,
                list) or isinstance(
                inputReynolds,
                np.ndarray):
            for re in inputReynolds:
                xf.find_coefficients(
                    airfoil=self.flapFoilName,
                    alpha=inputAlpha,
                    Reynolds=re,
                    iteration=inputIteration,
                    NACA=False,
                    delete=False,
                    dir=self.relFlapFoilDir() + '/')
        else:
            xf.find_coefficients(
                airfoil=self.flapFoilName,
                alpha=inputAlpha,
                Reynolds=inputReynolds,
                iteration=inputIteration,
                NACA=False,
                delete=False,
                dir=self.relFlapFoilDir())



def main():
    airf1 = rawFoil('EPA_0', 'G0', 'R0')
    print(airf1.flapFoilPath())
    airf1.loadFoil("FE_H5.dat")
    reArray = np.arange(1e5, 3e6, 1e5)
    alphaArray = np.arange(0, 30, 1)
    airf1.analyzeFoil(alphaArray, reArray, 20)

if __name__ == "__main__":
    main()
