##############################################################################
#  Unified RPA Primer & Molecular Beacon Designer for PfAgo / MB Diagnostics #
#  Miles Thompson, Pim Jitnavasathian & Madhu Prakash for Purdue iGEM 2020   #
#                                                                            #
#  Dependencies:                                                             #
#     Python 3.7.7                                                           #
#     selenium 3.141.0                                                       #
#     time, any version                                                      #
#     numpy, any version                                                     #
#     ChromeDriver: the latest stable version                                #
##############################################################################

# File Dependencies: 
#   A text file containing GenBank MT481992.1, assumed to be named Finalsequence.txt, in the working directory.
#   (Optional) A text file containing previous output of the calculateTargets function, to eliminate the need for excess webscraping in further runs. 

#############################################################################
# Code description: Input a genome, find optimal TtAgo guides, and for each #
#   Argonaute guide find the optimal primers for Recombinase Polymerase     #
#   Amplification (RPA) and optimal molecular beacons.                      #                
# Currently relies on user selection of a guide in the Executables section, #
#   followed by manual primer analysis in IDT's OligoAnalyzer. Manual       #
#   analysis of beacons on Quickfold server is also highly recommended.     #
#############################################################################

##############################
# SECTION 1: PACKAGE IMPORTS #
##############################

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time
import numpy as np

###############################
# SECTION 2: CLASS DEFINITION #
###############################

# Class that defines the mfoldSequence object. Each instance holds a sequence and other 
# information regarding a single molecular beacon as found on mfold. 
class mfoldSequence():
    # Sequence is the sequence of a given probe. Currently it's a list of chars. 
    # Hybridizations is an array of the indices of all hybridized bases in the sequence.
    # DG is a float representing the Gibbs Free Energy of the sequence.
    def __init__(self,sequence,hybridizations,DG):
        self.sequence = sequence
        self.hybridizations = hybridizations
        self.DG = DG

# Class that defines target sequences for primer design. Also houses all associated
# primers for an individual TtAgo cleavage system. 
class argonauteTarget():
    # Target is the specific portion of the sequence used to design Argonaute guides.
    # Extended target houses Target and has been extended equally in the 3' and 5' directions
    #   to enable RPA primer selection. 
    # gnGC measures the GC content of the gn guide, which is used as a metric for TtAgo target quality.
    # GT-GN are primers specific to the TtAgo cleavage of the target in dsDNA format. 
    # probe is the molecular beacon probe sequence that is housed within the target sequence. 
    # leftPrimer and rightPrimer are lists holding all primers found by primer3. Common primers share a common list index magnitude.
    def __init__(self,target,extendedTarget,GT,GF,GR,GN,probe,gnGC,leftPrimers=[],rightPrimers=[]):
        self.target = target
        self.extendedTarget = extendedTarget
        self.GT = GT
        self.GF = GF
        self.GR = GR
        self.GN = GN
        self.probe = probe
        self.gnGC = gnGC
        self.leftPrimers = leftPrimers
        self.rightPrimers = rightPrimers
        
# Unused class to hold a primer and its oligoAnalyzer output. 
# oA output retrieval currently not implemented. 
class rpaPrimer():
    def __init__(self,primer,output="NONE"):
        self.primer = primer
        self.output = output

####################################################
# SECTION 3: GLOBAL VARIABLE (& object) INITIATION #
####################################################

# Apparently this is poor programming practice. I do not care.
# Find here all parameters used on websites, as well as relevant file names.
# Change these variables to their desired values before running the code.

seqFileName = "Finalsequence.txt" # Name of the sequence file from which Argonaute guides are to be calculated. 
primerResultsFileName = "primerResults.txt" # Name of the file holding primer3 output
foldResultsFileName = "foldResults.txt" # Store .ct file outputs
beaconResultsFileName = "beaconResults.txt"
oaResultsFileName = "oligoResults.txt" # Name of the file holding oligoAnalyzer output 
finalResultsFileName = "primerFinalResults.txt"

# primer3 parameters
minPrimerLen = 30
maxPrimerLen = 36
optPrimerLen = 33
minPrimerGC = 30 # (I was told 20, but I think it was a typo.)
maxPrimerGC = 70
optPrimerGC = 50
minPrimerTM = 50
maxPrimerTM = 100
optPrimerTM = 75
prodRange = "100-200" # I don't know what the product range is or means. 
maxRepeat = 5 # Maximum tolerated mononucleotide repeat (inclusive?)

beaconNa = 0.015 # Hailey thinks 0.01
beaconMg = 0.005 # originally 0.01 but for some reason it was sitting at 0.005

# Length of beacon stems. If you wish to consider multiple lengths, you'll have 
# to change this value and rerun the code. 
stemLength = 4  # Default value: 5. Functional values: 4, 5. [6, 7 exceed quikfold's 1000-fold limit and would have to be broken up. That would be excessive use of a free service.]

# Temperature to fold the beacon at, in degrees C. 
temperature = 39  # Default value: 37. 

###################################
# SECTION 4: FUNCTION DEFINITIONS #
###################################

# Returns the target sequences for primer selection. fname is the name of the file housing the target genome.
def calculateTargets(fname):
    targets = []
    # I did not author this function. I'm not sure how it works, I just modified it for my purposes. It comes from Guide Program by Madhu Prakash and Pim Jitnavasathian.
    for nt in(16,17,18):# calculates guides of 16, 17, and 18nt. 
        infile = open(fname,'r')#opens user inputted file on read
        # initializing variables
        array = [] # houses targets
        extendedarray = [] # houses extended targets
        length = nt-16 # accounts for base pair shifts between 16,17,and 18nt. 
        for line in infile.readlines():# for everyline in the file
          for i in range(0,len(line)-(34+length)):#identifies 35,36, or 37nt sequences from the genom file from which gr,gn,gt,and gf will be designed from
              if line[i] == "t" and line[i+11] in ("g","c"): # checks for t at 5' of gt and a g/c at the 12th postion from the 5' of gt
                  if line[i+14+length] in ("g","c") and line[i+16 +length] == "t": # checks g/c at 12th position of gn from the 5' and t at the 5' of gf
                      if line[i+24+length] in ("g", "c") and line[i+25+length] == "a": #checks for g/ at 12th position of gr from 5' and a at the complement of the 5' of gn 
                          if line[i+27+length] in("g","c") and line[i+35+length] == "a": #checks for g/g at 12th position of gf from 5' and a at the complement of the 5' of gr
                              array.append(line[i:i+36+length]) #save resulting 35,36,or 37nt sequence in array
                              try: 
                                  extendedarray.append(line[i-82:i+36+length+82])
                              except: 
                                  extendedarray.append("ERROR: target too close to genome terminus.")
        for targetIndex in range(0,len(array)):
            gt = array[targetIndex][0:16+length]#identifies gt
            gf = array[targetIndex][16+length:32+(length*2)]#identifies gf
            gr_complement = array[targetIndex][20:36+length]#the "gr" in the array is the complement of the actual gr
            gn_complement = array[targetIndex][10:26+length]#the "gn" in the array is the complement of the actual gr
            g_content = gn_complement.count("g")#counts g's in gn complement
            c_content = gn_complement.count("c")#counts c's in gn complement
            gc_contentpercent = (g_content + c_content)/len(gn_complement)# calculates g/c percentage

            # identifies guide gr from gr complement
            gr1 = gr_complement.replace("t","A")
            gr2 = gr1.replace("g","C")
            gr3 = gr2.replace("a","T")
            gr_3_5 = gr3.replace("c","G") #gr from 3'--5'
            gr_5_3 = gr_3_5[::-1] #gr from 5'--3'

            # identifies guide gn from gn complement
            gn1 = gn_complement.replace("t","A")
            gn2 = gn1.replace("g","C")
            gn3 = gn2.replace("a","T")
            gn_3_5 = gn3.replace("c","G") #gn from 3'--5'
            gn_5_3 = gn_3_5[::-1] #gr from 5'--3'
            
            # checks if G/C content is less than .38 and greater than .12
            # populates finalarray with argonauteGuides objects
            if gc_contentpercent >= .12  and gc_contentpercent <= .38:
                # saves those guides along with MB probe and g/c percentage to a different array
                targets.append(argonauteTarget(array[targetIndex],extendedarray[targetIndex],gt.upper(),gf.upper(),gr_5_3.upper(),gn_3_5.upper(),gn_complement.upper(),gc_contentpercent))
        infile.close()
    return(targets) # return target array. 

# Void function selects primer sets from the online primer designer primer3 for each 
# target in the previously-calculated list of argonauteTargets. 
# Function will write all results to an outfile if given a name to write to,
#   and will read all results from a previously-generated file (bypassing the need for webscraping) 
#   if given an infile name. Do not pass both an outfile and infile name in the same run. 
def findPrimers(targetList, outFileName = "NONE", inFileName = "NONE"):
    if inFileName == "NONE": 
        driver = webdriver.Chrome() # Opens a Chrome window and instantiates a webdriver object to control it.
        driver.get("http://bioinfo.ut.ee/primer3/") # Navigate to the primer3 website.
        # Input all desired parameters. Done only once. 
        driver.find_element_by_name("PRIMER_MIN_SIZE").send_keys(Keys.BACKSPACE + Keys.BACKSPACE + str(minPrimerLen)) 
        time.sleep(0.25)
        driver.find_element_by_name("PRIMER_MAX_SIZE").send_keys(Keys.BACKSPACE + Keys.BACKSPACE + str(maxPrimerLen)) 
        time.sleep(0.25)
        driver.find_element_by_name("PRIMER_OPT_SIZE").send_keys(Keys.BACKSPACE + Keys.BACKSPACE + str(optPrimerLen))
        time.sleep(0.25)
        driver.find_element_by_name("PRIMER_MIN_TM").send_keys(Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + str(minPrimerTM))
        time.sleep(0.25)
        driver.find_element_by_name("PRIMER_MAX_TM").send_keys(Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + str(maxPrimerTM))
        time.sleep(0.25)
        driver.find_element_by_name("PRIMER_OPT_TM").send_keys(Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + str(optPrimerTM))
        time.sleep(0.25)
        driver.find_element_by_name("PRIMER_MIN_GC").send_keys(Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + str(minPrimerGC)) 
        time.sleep(0.25)
        driver.find_element_by_name("PRIMER_MAX_GC").send_keys(Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + str(maxPrimerGC))
        time.sleep(0.25)
        driver.find_element_by_name("PRIMER_OPT_GC_PERCENT").send_keys(Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + Keys.BACKSPACE + str(optPrimerGC))
        time.sleep(0.25)
        driver.find_element_by_name("PRIMER_MAX_POLY_X").send_keys(Keys.BACKSPACE + str(maxRepeat))
        time.sleep(0.25)
        driver.find_element_by_name("PRIMER_PRODUCT_SIZE_RANGE").clear()
        time.sleep(0.25)
        driver.find_element_by_name("PRIMER_PRODUCT_SIZE_RANGE").send_keys(prodRange)
        time.sleep(0.25)
        
        if outFileName != "NONE": # If you want to write primer3 output to a file
            outFileID = open(outFileName, "w")
            # Get results for each target in the list; use results to populate each target's primer lists. 
            for targetIndex in range(0, len(targetList)):
                results = ""
                sequenceBox = driver.find_element_by_name("SEQUENCE_TEMPLATE")
                sequenceBox.clear()
                sequenceBox.send_keys(targetList[targetIndex].extendedTarget) # find the site element used to enter the target sequence, then enter it. 
                time.sleep(0.25)
                # Pick some primers
                driver.find_element_by_name("Pick Primers").click() # Get our results. 
                results = str(driver.find_element_by_xpath("/html/body").text)
                outFileID.write(str(results))
                outFileID.write("END OF RESULTS FOR GIVEN PRIMER")
                results.split("/n")
                # Add a primer each time one is found. 
                for lineIndex in range(0,len(results)): 
                    line = results[lineIndex]
                    leftP = ""
                    rightP = ""
                    if "LEFT PRIMER" in line:
                        line = " ".join(line.split()).split(" ")
                        leftP = [line[-1]]
                        targetList[targetIndex].leftPrimers = targetList[targetIndex].leftPrimers + leftP # You CANNOT use the .append() method here. 
                    if "RIGHT PRIMER" in line:
                        line = " ".join(line.split()).split(" ")
                        rightP = [line[-1]]
                        targetList[targetIndex].rightPrimers = targetList[targetIndex].rightPrimers + rightP # You CANNOT use the .append() method here. 
                time.sleep(0.25)
                driver.back() # Apparently this method isn't the world's most reliable. Here I'm trying to backnavigate to the previous page, where all parameters are already input, to save runtime. 
            outFileID.close()
            driver.close()
        else: # If you do not want to write primer3 results to a file
            for targetIndex in range(0, len(targetList)):
                results = ""
                sequenceBox = driver.find_element_by_name("SEQUENCE_TEMPLATE")
                sequenceBox.clear()
                sequenceBox.send_keys(targetList[targetIndex].extendedTarget) # find the site element used to enter the target sequence, then enter it. 
                time.sleep(0.25)
                # Pick some primers
                driver.find_element_by_name("Pick Primers").click() # Get our results. 
                results = str(driver.find_element_by_xpath("/html/body").text).split("\n")
                # Add a primer each time one is found. 
                for lineIndex in range(0,len(results)): 
                    line = results[lineIndex]
                    leftP = ""
                    rightP = ""
                    if "LEFT PRIMER" in line:
                        line = " ".join(line.split()).split(" ")
                        leftP = [line[-1]]
                        targetList[targetIndex].leftPrimers = targetList[targetIndex].leftPrimers + rpaPrimer(leftP) # You CANNOT use the .append() method here. 
                    if "RIGHT PRIMER" in line:
                        line = " ".join(line.split()).split(" ")
                        rightP = [line[-1]]
                        targetList[targetIndex].rightPrimers = targetList[targetIndex].rightPrimers + rpaPrimer(rightP) # You CANNOT use the .append() method here. 
                time.sleep(0.25)
                driver.back() # Apparently this method isn't the world's most reliable. Here I'm trying to backnavigate to the previous page, where all parameters are already input, to save runtime. 
            driver.close()
    else: # If you want to use previously-generated primer results
        inFileID = open(inFileName, "r")
        fileContents = str(inFileID.read())
        allResults = fileContents.split("END OF RESULTS FOR GIVEN PRIMER") 
        for resIndex in range(0,len(allResults)):
            allResults[resIndex] = allResults[resIndex].split("\n")
        for targetIndex in range(0, len(targetList)):
            results = allResults[targetIndex]
            # Add a primer each time one is found. 
            for lineIndex in range(0,len(results)): 
                line = results[lineIndex]
                leftP = ""
                rightP = ""
                if "LEFT PRIMER" in line:
                    line = " ".join(line.split()).split(" ")
                    leftP = [line[-1]]
                    targetList[targetIndex].leftPrimers = targetList[targetIndex].leftPrimers + [rpaPrimer(leftP)] # You CANNOT use the .append() method here. 
                if "RIGHT PRIMER" in line:
                    line = " ".join(line.split()).split(" ")
                    rightP = [line[-1]]
                    targetList[targetIndex].rightPrimers = targetList[targetIndex].rightPrimers + [rpaPrimer(rightP)] # You CANNOT use the .append() method here.   
        inFileID.close()

# Function returns an array of all possible stems according to guidelines in
# "Molecular Beacon Sequence Design Algorithm", BioTechniques 34:68-73 (Jan. 2003)
def createStems(stemLength):
    arrayOfStems = []
    for letter1 in ['A','T','C','G']: 
        for letter2 in ['A','T','C','G']:
            for letter3 in ['A','T','C','G']:
                for letter4 in ['A','T','C','G']:
                    if stemLength > 4: 
                        for letter5 in ['A','T','C','G']: 
                            if stemLength > 5:
                                for letter6 in ['A','T','C','G']:
                                    if stemLength > 6: # 7 bp stem
                                        for letter7 in ['A','T','C','G']:    
                                            gcRatio = (int((letter1=='C' or letter1=='G'))+int((letter2=='C' or letter2=='G'))+int((letter3=='C' or letter3=='G'))+(letter4=='C' or letter4=='G')+int((letter5=='C' or letter5=='G'))+(letter6=='C' or letter6=='G')+(letter7=='C' or letter7=='G'))/7.0 
                                            if 0.7 <= gcRatio and 0.8 >= gcRatio:
                                                stem = letter1 + letter2 + letter3 + letter4 + letter5 + letter6 + letter7
                                                arrayOfStems.append(stem)
                                    else: # 6 bp stem
                                        letter7 = ""
                                        gcRatio = (int((letter1=='C' or letter1=='G'))+int((letter2=='C' or letter2=='G'))+int((letter3=='C' or letter3=='G'))+(letter4=='C' or letter4=='G')+int((letter5=='C' or letter5=='G'))+(letter6=='C' or letter6=='G'))/6.0    
                                        if 0.65 <= gcRatio and 0.85 >= gcRatio: # I loosened requirements here because otherwise a 6 bp stem would not be possible.
                                            stem = letter1 + letter2 + letter3 + letter4 + letter5 + letter6 + letter7
                                            arrayOfStems.append(stem)
                            else: # 5 bp stem
                                letter6 = ""
                                letter7 = ""
                                gcRatio = (int((letter1=='C' or letter1=='G'))+int((letter2=='C' or letter2=='G'))+int((letter3=='C' or letter3=='G'))+(letter4=='C' or letter4=='G')+int((letter5=='C' or letter5=='G')))/5.0     
                                if 0.7 <= gcRatio and 0.8 >= gcRatio:
                                    stem = letter1 + letter2 + letter3 + letter4 + letter5 + letter6 + letter7
                                    arrayOfStems.append(stem)
                    else: # 4 bp stem
                        letter5 = ""
                        letter6 = ""
                        letter7 = ""
                        gcRatio = (int((letter1=='C' or letter1=='G'))+int((letter2=='C' or letter2=='G'))+int((letter3=='C' or letter3=='G'))+(letter4=='C' or letter4=='G')+int((letter5=='C' or letter5=='G')))/4.0 # Determines fraction of nucleotides that are either g's or c's.
                        if 0.7 <= gcRatio and 0.8 >= gcRatio:
                            stem = letter1 + letter2 + letter3 + letter4 + letter5 + letter6 + letter7
                            arrayOfStems.append(stem)
    return arrayOfStems

# Returns array of the complements for all the stems. Credit to stack overflow.
def createAntiStems(arrayOfStems):
    revComplimentKey = arrayOfStems[0].maketrans("ACTG", "TGAC")
    antiStems = []
    for stem in arrayOfStems:
        antiStems.append(stem.translate(revComplimentKey)[::-1])
    return antiStems

# Populate allPossibleBeacons and allBeaconsList FOR A GIVEN PROBE. allPossibleBeacons will be input to Quikfold.
def generateBeacons(probeSequence,arrayOfStems,arrayOfAntiStems):
    quikfoldInput = ''
    allBeaconsList = []
    for index in range(0, len(arrayOfStems)):
        allBeaconsList.append(arrayOfStems[index] + probeSequence + arrayOfAntiStems[index])
        quikfoldInput = quikfoldInput + arrayOfStems[index] + probeSequence + arrayOfAntiStems[index] + ';'
    return allBeaconsList,quikfoldInput

# Determine if you need to fold. 
def checkIfAlreadyFolded():
    alreadyFolded = input("Do you already have a file for your beacon output? Input Yes or No. \n")
    if alreadyFolded == "Yes":
        return True
    if alreadyFolded == "No":
        return False

# Run all molecular beacon sequences through the Quikfold server. Potential to write them to a file so no duplicate folding will be necessary in the future.
def foldBeacons(quikfoldInput, outfileName = "NONE"):
    driver = webdriver.Chrome()
    # Go to the Quikfold server website where we'll fold all possible beacons.     
    driver.get("http://unafold.rna.albany.edu/?q=DINAMelt/Quickfold")
    multiSequenceBox = driver.find_element_by_name("seq") 
    multiSequenceBox.send_keys(quikfoldInput) # Enter sequence into box
    temperatureBox = driver.find_element_by_name('temp')
    temperatureBox.send_keys(Keys.BACKSPACE + Keys.BACKSPACE + str(temperature)) # Enter folding temperature 
    sodiumBox = driver.find_element_by_name('Sodium')
    sodiumBox.send_keys(Keys.BACKSPACE + str(beaconNa)) # Enter folding [Na+]
    magnesiumBox = driver.find_element_by_name('Magnesium')
    magnesiumBox.send_keys(Keys.BACKSPACE + str(beaconMg)) # Enter folding [Mg2+]
    driver.find_element_by_xpath('//*[@id="f"]/p[9]/input[1]').click() # Submit job
    driver.find_element_by_xpath('//*[@id="content"]/div/div[2]/div/div/p[4]/a[1]').click() # Get ct file
    foldOutput = driver.find_element_by_xpath('/html/body/pre').text
    if outfileName != "NONE": # If you want to add your .ct file results to a file for improved computational efficiency. 
        outfileID = open(outfileName, "w")
        outfileID.write(str(foldOutput))
        outfileID.close()
    return foldOutput 

# Void function outputs the fold results to a file in case you want to re-run without using Quikfold
def outputFolds(foldOutput,outfileName):
    outfileID = open(outfileName, "w")
    outfileID.write(str(foldOutput))
    outfileID.close()

# Splits output into a 3D array. First dimension represents individual beacons.
# Second dimension represents rows within the probe output.
# Third dimension represents the elements, by column, within each row. 
# dGs corresponding to each beacon are also returned. 
def splitString(string,lengthOfBeacon):
    stringArray = string.split(str(lengthOfBeacon) + " dG = ")
    dGs = [stringArray[index][0:stringArray[index].find(' ')] for index in range(1,len(stringArray))]
    dGs = [float(dGs[index]) for index in range(0,len(dGs))]
    stringArray = [stringArray[index].split("\n")[1:lengthOfBeacon+1] for index in range(1,len(stringArray))]
    for i in range(0,len(stringArray)):
        for j in range(0,len(stringArray[1])):
            stringArray[i][j] = stringArray[i][j].split(' ')
    return(stringArray, dGs)

# Given a specified beacon within foldOutput, return the a list of the bp 
#   indices at which it is hybridized. 
def getHybridizations(foldOutput,beaconIndex):
    indicesOfHybridizations = []
    for rowIndex in range(0,len(foldOutput[beaconIndex])):
        if int(foldOutput[beaconIndex][rowIndex][4]) != 0:
            indicesOfHybridizations.append(rowIndex)
    return(indicesOfHybridizations)

# Remove duplicate sequences; retain the minimum DG sequence if dupes are removed.
# This only works if the sequence list passed is sorted with minimum DGs first; this
# is the case in our instance (yay!). 
# APPARENTLY, THIS FUNCTION DOES NOT FUNCTION. WHICH IS NOT GOOD. 
def deDupe(sequences):
    deDuped = []
    for i in range(0, len(sequences)-1): 
        justSequences = [mfoldSequence.sequence for mfoldSequence in deDuped]
        if sequences[i].sequence not in justSequences: 
            deDuped.append(sequences[i]) 
    return(deDuped) 

# Codeblock to reformat the .ct file that Quikfold gave us into a list of 
# sequence objects, each representing a distinct [potential] beacon. 
def instantiateBeacons(foldOutput,lengthOfBeacon):
    # List to hold all mfold output sequence objects
    mfoldSequences = [] 
    dGs = []
    # Do some weird shenanigans that I don't remember to get nicely-formatted data.
    foldOutput, dGs = splitString(foldOutput,lengthOfBeacon)
    # Some more shenanigans, this time to put everything in an object. 
    for index in range(0, len(foldOutput)):
        mfoldSequences.append(mfoldSequence(list(np.transpose(foldOutput[index])[1]), getHybridizations(foldOutput, index), dGs[index]))    
    # Remove duplicate sequences, keeping only the most energetically favorable sequence. 
    mfoldSequences = deDupe(mfoldSequences) # this function was not individually tested yeet
    return mfoldSequences, dGs

# Takes the .ct output from Quikfold and determines which sequences are acceptable molecular beacons. 
# Does this by making sure that only the stem nucleotides are hybridized. 
# POTENTIALLY LOGICALLY ERRONEOUS??? 
def selectForGoodBeacons(mfoldSequences, lengthOfBeacon):
    goodBeaconIndices = []
    for beaconIndex in range(0,len(mfoldSequences)):
        indicesOfHybridizations = mfoldSequences[beaconIndex].hybridizations
        correctNumberOfHybridizations = (len(indicesOfHybridizations) == 2 * stemLength)
        correctIndicesOfHybridizations = (indicesOfHybridizations == list(range(stemLength)) + list(range(lengthOfBeacon-stemLength, lengthOfBeacon))) # um what the falcon is this
        if correctNumberOfHybridizations and correctIndicesOfHybridizations: 
            # Make sure your logic is correct about the output counting a hybridization for each member of hybridized pairs.
            goodBeaconIndices.append(beaconIndex)
    return(goodBeaconIndices)

# Get rid of non alphanumerics in my output. Credit to DrAl of Stack Overflow.
def ExtractAlphanumeric(InputString):
    from string import ascii_letters, digits
    return "".join([ch for ch in InputString if ch in (ascii_letters + digits)])

def outputPrimerResults(targetList, outFileName, showAll = False):
    outFileID = open(outFileName, "w")
    outFileID.write("Results of Primer Design:\n")
    if showAll == False:
        outFileID.write("(This only displays the optimal primers according to primer3)\n")
    else:
        outFileID.write("\n")
    for argTarg in targetList:
        outFileID.write("\nTarget: " + argTarg.target + "\n")
        if showAll == True:
            for primerIndex in range(0,len(argTarg.leftPrimers)):
                outFileID.write("Left Primer " + str(primerIndex+1) + ": " + str(argTarg.leftPrimers[primerIndex].primer) + ";\n")
                outFileID.write("Right Primer " + str(primerIndex+1) + ": " + str(argTarg.rightPrimers[primerIndex].primer) + "\n")
        else:
            for primerIndex in range(0,1):
                outFileID.write("Left Primer " + str(primerIndex+1) + ": " + "".join(argTarg.leftPrimers[primerIndex].primer) + "\n")
                outFileID.write("Right Primer " + str(primerIndex+1) + ": " + "".join(argTarg.rightPrimers[primerIndex].primer) + "\n")          
    outFileID.close()
    
def outputBeaconResults(bestSequences, outFileName, numToShow, probeSequence):
    outFileID = open(outFileName, "w")
    outFileID.write("For the probe sequence " + str(probeSequence) + ", the " + str(numToShow) + " best beacons, in descending order, are:\n")
    for beaconIndex in range(0,min(numToShow,len(bestSequences))):
        outFileID.write("\n" + bestSequences[beaconIndex] + "\n")
    print("\nPlease find beacon selection results in the following file: " + str(outFileName) + "\nMore detail is held in the goodBeacons variable.")
    outFileID.close()

##############################################################################
##############################################################################
########################## SECTION 5: EXECUTABLES ############################
##############################################################################
##############################################################################

targets = calculateTargets(seqFileName) # Find RPA amplification targets that will be cleaved by Argonaute. This function is adapted from Pim's and Madhu's code and I can't speak to the biology underpinning it. 
targets.sort(key = lambda argonauteTarget: argonauteTarget.gnGC) # Sort targets by increasing gc percents of their gn sequences. 
findPrimers(targets, inFileName = primerResultsFileName) # Find optimal primers for each amplification target using primer3. 
arrayOfStems = createStems(stemLength) # Generate all possible molecular beacon stems of the given length. Generates the ssDNA portion of the stem that attaches to the 5' end of the probe.
arrayOfAntiStems = createAntiStems(arrayOfStems) # Generate the complement to the previous stems, which will attach at the 3' end of the probe.
# IMPORTANT: User input required within the program! 
probeSequence = targets[5].probe # This probe is the only one that corresponds to both good RPA primers and a functional molecular beacon. You will manually select one of the guides for each run of the program and inspect the generated primers and beacon(s).
lengthOfBeacon = len(probeSequence) + 2 * stemLength 
alreadyFolded = checkIfAlreadyFolded() # Have you already saved the beacon folding results to a file? If so, it will read from the file. Else it will webscrape from QuickFold. 
foldOutput = ""
if not alreadyFolded:
    allBeaconsList, quikfoldInput = generateBeacons(probeSequence,arrayOfStems,arrayOfAntiStems) # Make a list of all beacons (left stem + probe + right stem) 
    foldOutput = foldBeacons(quikfoldInput, outfileName = foldResultsFileName) # Webscrape from Quickfold to find good molecular beacons. 
else: 
    foldResultsID = open(foldResultsFileName, "r") # Read from file with webscraping results. 
    foldOutput = foldResultsID.read() 
    foldResultsID.close() 
mfoldSequences, dGs = instantiateBeacons(foldOutput,lengthOfBeacon) # Create a beacon object for each beacon sequence by taking into account folding information. 
goodBeaconIndices = selectForGoodBeacons(mfoldSequences,lengthOfBeacon) # Find the indices of good molecular beacons.
goodBeacons = [mfoldSequences[i] for i in goodBeaconIndices] # Generate a list of the good molecular beacon sequences. 
for beacon in goodBeacons:
   beacon.sequence = ExtractAlphanumeric(beacon.sequence) # Formatting
goodBeacons.sort(key = lambda mfoldSequence: mfoldSequence.DG) # Sort by decreasing spontaneity of the beacon's most spontaneous structure
bestSequences = [ExtractAlphanumeric(beacon.sequence) for beacon in goodBeacons] # Sorted optimal sequences in string format
try:
    if goodBeaconIndices[0] == goodBeaconIndices[1]:
        print("\nNo good beacons found.")
    else:
        outputBeaconResults(bestSequences,beaconResultsFileName,5,probeSequence)
except:
    print("No good beacons found for that probe. Sorry.") # If you get this output you cannot use the Argonaute / RPA target you selected. Re-run. 