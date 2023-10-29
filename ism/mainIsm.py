
# MAIN FUNCTION TO CALL THE ISM MODULE

from ism.src.ism import ism

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = r'C:\Users\diego\PycharmProjects\pythonProjectNewGIT\auxiliary'
indir = r"C:\\Users\\diego\\PycharmProjects\\EODP_TER_2021\\EODP-TS-ISM\\input\\gradient_alt100_act150" # small scene
# indir = r"C:\\Users\\diego\\PycharmProjects\\EODP_TER_2021\\EODP-TS-E2E\\sgm_out"
outdir = r"C:\\Users\\diego\\PycharmProjects\\EODP_TER_2021\\EODP-TS-ISM\\myoutput\\TestDetectionModule"
# outdir = r"C:\\Users\\diego\\PycharmProjects\\EODP_TER_2021\\EODP-TS-ISM\\myoutput\\newOutput"





# Initialise the ISM
myIsm = ism(auxdir, indir, outdir)
myIsm.processModule()
