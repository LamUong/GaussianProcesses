import sys, os, subprocess
import matplotlib.pyplot as plt
from matplotlib import gridspec
cwd = os.getcwd()
sys.path.append(cwd + '/../tip3p/')
import datagenerator

# ----- globals -----
moleculeFName = "molecule.init"
Action = False

def processDatapoint(datapoint, mainDir, inputFName, fileIDdata, fileIDxyz):
   
    global Action

    params = datapoint[0]
    coords = datapoint[1]
    # convert to angstroms
    for item in range(len(coords)):
        for jtem in range(len(coords[item])):
            coords[item][jtem] /= datagenerator.ang2bohr

    labelStr = getParamLabel(params)

    # define and create job directory
    jobDir = mainDir + "/" + labelStr

    # submission script name
    qsubFile = inputFName + ".qsub"

    # the goal is to perform all operations that are necessary to get the final energy
    # all encountered problems should be resolved along the way automatically
    if not os.path.exists(jobDir):

        # create a new directory and submit job
        os.makedirs(jobDir)

        # create link to the input file
        fname = jobDir + "/" + inputFName
        templateInp =  "../" + inputFName
        os.symlink(templateInp, fname)

        # create link to the submission file
        fname = jobDir + "/" + qsubFile
        templateInp =  "../" + qsubFile
        os.symlink(templateInp, fname)

        # create file with molecule description
        fname = jobDir + "/" + moleculeFName 
        fileID = open( fname, 'w' )
        coordStr = getCoordStr(0, coords)
        fileID.write(coordStr)
        fileID.close()

        coordStr = getCoordStr(2, coords)
        coordStr = fileIDxyz.write(coordStr)

        submitJob(Action, jobDir, qsubFile)

    else: # job directory exists

        # check if the energy is already calculated
        # if so collect, if not re-submit
        print ("Directory %s exists!" % jobDir)

        outfile = inputFName + ".out"
        if not os.path.exists(outfile):
            # (re-)submit
            submitJob(Action, jobDir, qsubFile)
        else:
            # analyze output
            energy = analyzeJob(Action, jobDir, outfile)
            if (energy < -0.00001):
                # energy is OK, make a record
                coordStr = getCoordStr(1,coords)
                coordStr = fileIDdata.write(coordStr)
                coordStr = getCoordStr(2,coords)
                coordStr = fileIDxyz.write(coordStr)
            else:
                # re-submit
                submitJob(Action, jobDir, qsubFile)

def analyzeJob(action, jobDir, outfile):

    if (action):
        command = "cd %s; grep 'CCSD(T) Total Energy' %s | awk '{print $5}' " % (jobDir, outfile)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        proc_stdout = process.communicate()[0].strip()
        energy = float(proc_stdout)
    else:
        print ( "Get E: cd %s; grep 'CCSD(T) Total Energy' %s | awk '{print $5}' " % (jobDir, outfile) )
        energy = 0.0

    return energy

def submitJob(action, jobDir, qsubfile):

    if (action):
        command = "cd %s; qsub %s" % (jobDir, qsubfile)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        proc_stdout = process.communicate()[0].strip()
        print (proc_stdout)
    else:
        print ( "Submit: cd %s; qsub %s" % (jobDir, qsubfile) )

def getCoordStr(formattype, coords):
    #formattype: 0 - qchem
    #            1 - fitting database
    #            2 - xyz

    natoms = len(coords)
    if ( natoms%datagenerator.atoms_per_mol != 0 ):
        print ("Corrupt data")
        exit(2)
                      
    if (formattype==0):
        coordStr = '$molecule\n0 1\n'
    elif (formattype==2):
        coordStr = "%d\ncomment\n" % natoms 

    for iatom in range(0, natoms):
        #if ( iatom%datagenerator.atoms_per_mol == 0 ):
        #    coordStr += '--\n0 1\n'
        atomStr = atomSymbol(iatom)
        if (formattype==1):
            atomStr += ":"
        for icoord in range(3):
            atomStr += "%20.10f" % coords[iatom][icoord]
        atomStr += "\n"
        coordStr += atomStr

    if (formattype==0):
        coordStr += '$end\n'

    return coordStr

def atomSymbol(iatom):
    if ( iatom%datagenerator.atoms_per_mol !=0 ):
        symbol = "%3s" % "H"
    else:
        symbol = "%3s" % "O"
    return symbol

def getParamLabel(params):

    nparams = len(params)
    labelStr = "Q"
    for iparam in range(nparams):
        if (iparam == 0 or iparam == 3 or iparam == 4):
            # these are distances in bohr
            value = params[iparam]/datagenerator.ang2bohr
            iparamStr = "-%07.3f" % value
        else:
            # these are angles in radians
            iparamStr = "-%07.3f" % params[iparam]
        labelStr += iparamStr

    return labelStr

def main(): 

    script, infile = sys.argv

    print ("Input file %s" % infile)
    if not os.path.exists(infile):
        print ("File %s does not exists!" % infile)
        exit(2)

    qsubfile = infile + ".qsub"
    if not os.path.exists(qsubfile):
        print ("File %s does not exists!" % qsubfile)
        exit(2)

    datafile = infile + ".data"
    if os.path.exists(datafile):
        print ("Error: File %s exists!" % datafile)
        exit(2)

    xyzfile = infile + ".xyz"
    if os.path.exists(xyzfile):
        print ("Error: File %s exists!" % xyzfile)
        exit(2)

    # create files to write results
    fileIDdata = open( datafile, 'w' )
    fileIDxyz = open( xyzfile, 'w' )

    QCinpFName = os.path.basename(infile)
    mainDir = os.path.dirname(infile)

    print("Generating data")
    data = datagenerator.main(energyRequired=False)

    npoints = len(data)
    for point in range(0,npoints):
        processDatapoint(data[point], mainDir, QCinpFName, fileIDdata, fileIDxyz)

    fileIDdata.close()
    fileIDxyz.close()

main()

