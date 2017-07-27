import os
from multiprocessing import Pool
import slurm

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"

def processHelper(fileName):
    options = "-WAIT -fillsidechains -f 3 -fix -samplewater -delwater_hbond_cutoff 2 -keepfarwat -captermini"  
    command = "{}/utilities/prepwizard {} -j temp-{} ../stripped/{}.mae {}.mae".format(SCHRODINGER, options, fileName, fileName, fileName)
    slurm.salloc(command, "1", "5:00:00") #(jobString, numProcessors, timeLimit)
    
    return fileName

def processSuccess(structure):
    fileName = structure+".mae"
    return fileName in os.listdir(".") #We're currently in processed

def process():   
    toProcess = [os.path.splitext(f)[0] for f in os.listdir("stripped") if not os.path.isfile('processed/'+f)]
    print toProcess
    os.system('mkdir -p processed')
    os.chdir('processed')
    #Each thread has a very lightweight job (waiting on the salloc signal) so just launch len(currentFiles) threads
    pool = Pool(len(toProcess))
    
    #Keep submitting jobs until they're all successful
    while len(toProcess) != 0:
        currentlyProcessing = toProcess
        toProcess = []

        for finishedStruct in pool.imap_unordered(processHelper, currentlyProcessing):
            if not processSuccess(finishedStruct):
                toProcess.append(finishedStruct)
                print("{} did not process successfully - resubmitting".format(finishedStruct))

    os.system('rm temp* rm *missing* rm *.out')
    os.chdir('..')