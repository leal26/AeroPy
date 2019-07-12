from rapidboom import AxieBump
from weather.boom import read_input
from weather.scraper.twister import process_data
import platform


def calculate_loudness(bump_function):
    # Bump design variables
    location = 12.5

    # Flight conditions inputs
    alt_ft = 50000.

    # Setting up for
    CASE_DIR = "./"  # axie bump case
    PANAIR_EXE = 'panair.exe'
    SBOOM_EXE = 'sboom_windows.dat.allow'

    # Run
    # axiebump = AxieBump(CASE_DIR, PANAIR_EXE, SBOOM_EXE) # for standard atmosphere
    axiebump = AxieBump(CASE_DIR, PANAIR_EXE, SBOOM_EXE, altitude=alt_ft,
                        deformation='custom')
    axiebump.MESH_COARSEN_TOL = 0.00045
    axiebump.N_TANGENTIAL = 20
    loudness = axiebump.run(bump_function)

    return loudness
