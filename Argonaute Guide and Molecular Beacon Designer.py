######################################################
# Argonaute Guide & Molecular Beacon Designer        #
# Miles Thompson, Pim Jitnavasathian & Madhu Prakash #
# Purdue iGEM Summer Team 2020                       #
######################################################

# Input: A DNA sequence for use in an Argonaute / molecular beacon detection 
#   system. This code was designed to use GenBank: MT481992.1, but should be 
#   robust. 

# Output: The optimal Argonaute guide sequences, as well as the optimal 
#   molecular beacon sequence, for this detection system. 

# Notable Editable Parameters: See section 4. 

# There are separate files that can independently optimize for (1) Ago guide
#   sequences given a genome and (2) molecular beacons given a probe sequence.
#   They are available on request, and at the time of writing they were called 
#   Guide Program.py and MolecularBeaconDesignerV3.py, respectively. 

####################
# SECTION 1: SETUP #
####################

# Setup instructions for the (required) selenium package are available here:
# https://selenium-python.readthedocs.io/installation.html

# For Windows machines, here are the basics: 
# 1.) Get python running in your command terminal. (Put your python instance
#     in your cmd path. This means you need your "python" program in C:\\users\\yourname),
#     or you must be able to navigate to python using the cd command.
# 2.) Use the command "pip install selenium" to install the selenium package. 
#     If this doesn't work and you have Anaconda, try "conda install selenium".
# 3.) Download the most recent non-beta Chrome webdriver for your system. It
#     must be saved to the same folder as your "python" program. Downloads are 
#     available here: 
#     https://chromedriver.chromium.org/

##############################
# SECTION 2: PACKAGE IMPORTS #
##############################

import selenium # webdriver package to access primer3, quickfold, & other sites
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import random # enables us to be sneaky
import time # also enables sneakiness
import numpy as np 

###############################
# SECTION 3: CLASS DEFINITION #
###############################

# Class that defines the mfoldSequence object. This will to make our mfold 
# output nicer. It also makes everything a little cleaner. 
class mfoldSequence():
    # Sequence is the sequence of a given probe. Currently it's a list of chars. 
    # Hybridizations is an array of the indices of all hybridized bases in the sequence.
    # DG is a float representing the Gibbs Free Energy of the sequence.
    def __init__(self,sequence,hybridizations,DG):
        self.sequence = sequence
        self.hybridizations = hybridizations
        self.DG = DG

# Class that defines the argonauteGuide object. It holds a GT (5'-3'),
# GF (5'-3'), GR (5'-3'), and GN (3'-5') sequence, as well as a molecular 
# beacon probe. The GC% of the GN / probe is also stored.
class argonauteGuides():
    def __init__(self,GT,GF,GR,GN,probe,gcPercent):
        self.GT = GT
        self.GF = GF
        self.GR = GR
        self.GN = GN
        self.probe = probe
        self.gcPercent = gcPercent
        
####################################################
# SECTION 4: GLOBAL VARIABLE (& object) INITIATION #
####################################################

# Apparently this is poor programming practice. I do not care.
# Change these variables to their desired values before running the code.

# Length of beacon stems. If you wish to consider multiple lengths, you'll have 
# to change this value and rerun the code. 
stemLength = 5  # Default value: 5. Functional values: 4, 5, [6?, 7? (exceed quikfold's 1000-fold limit)]

# Temperature to fold the beacon at, in degrees C. 
temperature = 37  # Default value: 37. Functional values: figure it out. 

driver = webdriver.Chrome() # Instantiate a chrome webdriver object to control the web. 

# Desired [Na+] to fold the beacon at. Modify here, not in user input.
naConc = 0.01 # Default value: 0.01. Functional values: 0.01-?

# Desired [Mg2+] to fold the beacon at. Modify here, not in user input.
mgConc = 0.001 # Default value: 0.001. Functional values: ?-?

###################################
# SECTION 5: FUNCTION DEFINITIONS #
###################################

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
    return(arrayOfStems)

# Returns array of the complements for all the stems. Credit to stack overflow.
def createAntiStems(arrayOfStems):
    revComplimentKey = arrayOfStems[0].maketrans("ACTG", "TGAC")
    antiStems = []
    for stem in arrayOfStems:
        antiStems.append(stem.translate(revComplimentKey)[::-1])
    return antiStems

# Splits output into a 3D array. First dimension represents individual beacons.
# Second dimension represents rows within the probe output.
# Third dimension represents the elements, by column, within each row. 
# dGs corresponding to each beacon are also returned. 
def splitString(string):
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

# Takes the .ct output from Quikfold and determines which sequences are acceptable molecular beacons. 
# Does this by making sure that only the stem nucleotides are hybridized. 
# POTENTIALLY LOGICALLY ERRONEOUS??? 
def selectForGoodBeacons(mfoldSequences):
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

# function to return the complement of a DNA sequence. Not currently used.
def complement(sequence):
    sequence = sequence.upper()
    complementKey = sequence.maketrans("ACTG", "TGAC")
    comp = ""
    comp = comp + sequence.translate(complementKey)
    return comp

# function to return the reverse complement of a DNA sequence. Not currently used.
def reverseComplement(sequence):
    revComplimentKey = sequence.maketrans("ACTG", "TGAC")
    revcomp = ""
    revcomp = revcomp + sequence.translate(revComplimentKey)[::-1]
    return revcomp
        
## Function to calculate guides##
def Calculateguide():
    guides = []
    for nt in(16,17,18):# calculates guides of 16, 17, and 18nt. 
        infile = open(fname,'r')#opens user inputted file on read
        # initializing variables
        array = [] 
        finalarray = [] # will hold all guide sequences
        length = nt-16 # accounts for base pair shifts between 16,17,and 18nt. 
        for line in infile.readlines():# for everyline in the file
          for i in range(0,len(line)-(34+length)):#identifies 35,36, or 37nt sequences from the genom file from which gr,gn,gt,and gf will be designed from
              if line[i] == "t" and line[i+11] in ("g","c"): # checks for t at 5' of gt and a g/c at the 12th postion from the 5' of gt
                  if line[i+14+length] in ("g","c") and line[i+16 +length] == "t": # checks g/c at 12th position of gn from the 5' and t at the 5' of gf
                      if line[i+24+length] in ("g", "c") and line[i+25+length] == "a": #checks for g/ at 12th position of gr from 5' and a at the complement of the 5' of gn 
                          if line[i+27+length] in("g","c") and line[i+35+length] == "a": #checks for g/g at 12th position of gf from 5' and a at the complement of the 5' of gr
                              array.append(line[i:i+36+length]) #save resulting 35,36,or 37nt sequence in array
        # sorts through the arrays and respectively finds gr,gt,gn,and gf
        for x in array:
            gt = x[0:16+length]#identifies gt
            gf = x[16+length:32+(length*2)]#identifies gf
            gr_complement = x[20:36+length]#the "gr" in the array is the complement of the actual gr
            gn_complement = x[10:26+length]#the "gn" in the array is the complement of the actual gr
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
                guides.append(argonauteGuides(gt.upper(),gf.upper(),gr_5_3.upper(),gn_3_5.upper(), gn_complement.upper(), gc_contentpercent))
        infile.close()
    guides.sort(key = lambda argonauteGuides: argonauteGuides.gcPercent)
    return(guides) # return guides array sorted by increasing gc percents.)

#########################
# SECTION 6: USER INPUT #
#########################
    
# No, Bill Crum. I will not make a function for each user input. 

# Sequence to enter into primer3 / for use as probe. 
fname = input("Enter the file name of the sequence / genome you wish to detect with an Ago / MB system. Ensure that you add the .txt at the end:\n" ) #file must be lowercase no spaces and should be 1 big string.
guideSequences = Calculateguide()

prospectiveSequence = guideSequences[0].probe

useGivenSequence = input("I've found that " + str(prospectiveSequence) + " is likely the optimal probe sequence. If you want to use it, enter 1. If you want to input your own probe sequence, enter 0.\n")

if not useGivenSequence:
    probeSequence = input("Alright, enter your probe sequence. No spaces. \n").upper()

else:
    probeSequence = str(prospectiveSequence)
    
print("\nJust a moment, I'm running hundreds of sequences through online folding software.")
    
##########################
# SECTION 7: EXECUTABLES #
##########################
    
#/////////////////////////////////////////////////////////////////////////////#
# Stem-finding codeblock. Find stems. Yay.

arrayOfStems = createStems(stemLength)
arrayOfAntiStems = createAntiStems(arrayOfStems)

#/////////////////////////////////////////////////////////////////////////////#
# This codeblock formats the necessary input for the DNA folding software.

# String containing all possible beacon sequences.
allPossibleBeacons = ''
allBeaconsList = []

# Populate allPossibleBeacons; separate each sequence with a semicolon.
# I do this because it's the required input format of QuikFold. 
for index in range(0, len(arrayOfStems)):
    allBeaconsList.append(arrayOfStems[index] + probeSequence + arrayOfAntiStems[index])
    allPossibleBeacons = allPossibleBeacons + arrayOfStems[index] + probeSequence + arrayOfAntiStems[index] + ';'
    
lengthOfBeacon = len(probeSequence) + 2 * stemLength # I mean... just think about it.

#/////////////////////////////////////////////////////////////////////////////#
# Here we interact with the DNA folding software Quikfold. It outputs a .ct 
# file which we will interpret later on. 

# Go to the Quikfold server website where we'll fold all possible beacons.     
driver.get("http://unafold.rna.albany.edu/?q=DINAMelt/Quickfold")
multiSequenceBox = driver.find_element_by_name("seq") 
multiSequenceBox.send_keys(allPossibleBeacons) # Enter sequence into box
temperatureBox = driver.find_element_by_name('temp')
temperatureBox.send_keys(Keys.BACKSPACE + Keys.BACKSPACE + str(temperature)) # Enter folding temperature 
sodiumBox = driver.find_element_by_name('Sodium')
sodiumBox.send_keys(Keys.BACKSPACE + str(naConc)) # Enter folding [Na+]
magnesiumBox = driver.find_element_by_name('Magnesium')
magnesiumBox.send_keys(Keys.BACKSPACE + str(mgConc)) # Enter folding [Mg2+]
driver.find_element_by_xpath('//*[@id="f"]/p[9]/input[1]').click() # Submit job
driver.find_element_by_xpath('//*[@id="content"]/div/div[2]/div/div/p[4]/a[1]').click() # Get ct file

# Output of the DNA folding software QuikFold
foldOutput = driver.find_element_by_xpath('/html/body/pre').text

#/////////////////////////////////////////////////////////////////////////////#
# Codeblock to reformat the .ct file that Quikfold gave us into a list of 
# sequence objects, each representing a distinct [potential] beacon. 

# List to hold all mfold output sequence objects
mfoldSequences = [] 

# Do some weird shenanigans that I don't remember to get nicely-formatted data.
foldOutput, dGs = splitString(foldOutput)

# Some more shenanigans, this time to put everything in an object. 
for index in range(0, len(foldOutput)):
    mfoldSequences.append(mfoldSequence(list(np.transpose(foldOutput[index])[1]), getHybridizations(foldOutput, index), dGs[index]))

# Remove duplicate sequences, keeping only the most energetically favorable sequence. 
mfoldSequences = deDupe(mfoldSequences) # this function was not individually tested yeet

#/////////////////////////////////////////////////////////////////////////////#
# Select for all sequences that are potentially-functional beacons.

# Indices of foldOutput which correspond to sequences within foldOutput with
# the required folding structure. 
goodBeaconIndices = selectForGoodBeacons(mfoldSequences) 

#/////////////////////////////////////////////////////////////////////////////#
# Make a list of the best beacons sorted by minimum DG.

goodBeacons = [mfoldSequences[i] for i in goodBeaconIndices]

goodBeacons.sort(key = lambda mfoldSequence: mfoldSequence.DG)

bestSequences = [ExtractAlphanumeric(beacon.sequence) for beacon in goodBeacons]

"""
bestBeaconIndices = []
# Select the minimum 10 dGs from the good beacons.
dGs_copy = dGs 
for index in range(0, 10):
    for dGfixer in range(0,len(dGs_copy)):  # only consider the dGs of good beacons in optimization
        try:
            inTheList = goodBeaconIndices.index(dGfixer)
        except ValueError:
            dGs[dGfixer] = 100          
    bestBeaconIndices.append(dGs_copy.index(min(dGs_copy))) # add the best beacon first
    dGs_copy[dGs_copy.index(min(dGs_copy))] = 100 # remove the best beacon from consideration

# hey kid want the best sequence of course you do well here it is
bestSequences = [ExtractAlphanumeric(mfoldSequences[beaconIndex].sequence) for beaconIndex in bestBeaconIndices]
bestSequence = bestSequences[0]
"""

#############################
# SECTION 8: OUTPUT RESULTS #
#############################

if goodBeaconIndices[0] == goodBeaconIndices[1]:
    print("\nNo good beacons found.")
    
else:
    print("\nPlease run my top results through Quikfold and see how you like them. Here are the 10 best beacons, ordered by increasing DG:\n")
    print(bestSequences[0:9])
    print("\nSome tips: Avoid sequences whose most spontaneous structures have very close DGs. Preferably, we would like a single structure with a low DG. I may adjust the code to take this into account later.")
    
satisfied = input("\nWere you completely satisfied with one of these beacons? Enter 1 for yes and 0 for no.")

if not satisfied:
    print("\nPlease select a different probe sequence from this list and rerun the code. If you've done them all, try changing the stem length.")
    print([guideSequences[i].probe for i in range(0,len(guideSequences))])


else:
    satisfactorySequence = input("\nGreat! Please let me know which one:")
    print("\nHere's all the info you might need:")
    for index in range(0,len(guideSequences)): 
        if guideSequences[index].probe == satisfactorySequence:
            print("GT: " + guideSequences[index].GT + "\nGF: " + guideSequences[index].GF + "\nGR: " + guideSequences[index].GR + "\nGN: " + guideSequences[index].GN + "\nProbe: " + guideSequences[index].probe + "\nGC percent: " + str(guideSequences[index].gcPercent))

# TO-DO:
# ENSURE THAT THE PROBE SEQUENCES ARE ACTUALLY PROBE SEQUENCES BY CROSS REFERENCING WITH ORIGINAL GUIDE PROGRAM [strictly necessary!!!]
# PRIORITIZE SEQUENCES WITH NO DUPLICATES [hard]

###########################
# SECTION 9: HOUSEKEEPING #
###########################

driver.close() # Close the Chrome window and terminate the driver. 
driver.quit()  # Probably overkill but I don't like background processes. 

# HEY YOU. Read this. 
# If there is a runtime error, run the above two methods in the console. 