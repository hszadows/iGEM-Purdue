##################################
# Molecular Beacon Designer V2.0 #
# Miles Thompson                 #
# Purdue iGEM Summer Team 2020   #
##################################

# Input: DNA sequence to screen for the optimal subsequence for use as 
#   a molecular beacon probe.
# Output: I'm not quite sure what all it'll output. Certainly some beacon 
#   candidates at the very least. Currently it outputs what it thinks is the
#   best molecular beacon. bestBeaconIndices holds the indices of the best 
#   beacons in descending order of goodness. Index through mfoldSequences to 
#   check out the possible beacons and screen the output for the beacon you 
#   prefer.
# It is fairly common for the code not to find a good beacon. If this is the 
#   case, bestBeaconIndices will be populated exclusively by 0s and the code 
#   will output the first candidate sequence from mfoldSequences. This does 
#   not mean you should use this as a beacon. This means you should get a 
#   different internal oligo. 

# EVERY DNA SEQUENCE IS 5' TO 3' UNLESS OTHERWISE STATED.

# A note on running this code: I seriously doubt your python compiler / IDE 
# has access to selenium, which is required for my code. You'll need to 
# enable python to run in your Windows command terminal (look it up), then use
# the command "pip install selenium". If this doesn't work and you have 
# Anaconda (Spyder), you could try "conda install selenium". 

import selenium # webdriver package to access primer3, quickfold, & other sites
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import random
import time
import numpy as np

# The code currently works only for 4 and 5 bp stems. 6 bp stems don't work 
#   based on my beacon quality constraints. 7 bp stems make runtime excessive.
#   I could expand to 6 by loosening my constraints a little. 
# The code ONLY considers beacons with stems of ONE length at a time. If you 
#   want to test 4 and 5 bp stems, you'll have to run twice. 
stemLength = 5

# Temperature to fold the beacon at. In degrees C. 
temperature = 64

driver = webdriver.Chrome() # Open a Chrome window to control the web

# Sequence to enter into primer3. Will give us a probe sequence.
initialSequence = input("Enter a DNA sequence. No spaces, capitals only. \n") 
probeSequence = ""

dontNeedPrimer3 = input("Input 1 if this is the desired probe sequence. Else input 0. \n")

# If you need to use primer3 to find a probe sequence this code block does so.
if(not bool(int(dontNeedPrimer3))): 
    # Open primer3 in Chrome
    driver.get("http://bioinfo.ut.ee/primer3/")
    time.sleep(random.random()) # Wait for up to 1 second between commands. This helps 
                         # prevent the site from thinking we're just a bot and
                         # getting mad.
    
    # Select desired site settings (select for probe, not left or right primers)
    # These are specifically clicking checkboxes. 
    # NOTE: If you want to use specific temperatures / other conditions you'll need to do that here.
    driver.find_element_by_name("MUST_XLATE_PRIMER_PICK_INTERNAL_OLIGO").click()
    time.sleep(random.random())
    driver.find_element_by_name("MUST_XLATE_PRIMER_PICK_LEFT_PRIMER").click()
    time.sleep(random.random())
    driver.find_element_by_name("MUST_XLATE_PRIMER_PICK_RIGHT_PRIMER").click()
    
    # Enter the sequence.
    sequenceBox = driver.find_element_by_name("SEQUENCE_TEMPLATE")
    sequenceBox.send_keys(initialSequence)
    
    # Get probe sequence optimization results.
    driver.find_element_by_name("Pick Primers").click()
    
    # I am lazy. 
    probeSequence = input("Please do me a solid and paste the first listed internal oligo sequence here: \n")

else: # If you wanted to just use the originally input sequence as your probe 
    probeSequence = initialSequence
 
#
    #
#       #
    #       #
# FIND STEMS HERE
    #       #
#       #
    #
#
    
# Function returns an array of all possible stems according to guidelines in
# "Molecular Beacon Sequence Design Algorithm", BioTechniques 34:68-73 (Jan. 2003)
def createStems(stemLength):
    arrayOfStems = []
    for letter1 in ['A','T','C','G']: # Remove last G to prevent 5' Guanine? 
        for letter2 in ['A','T','C','G']:
            for letter3 in ['A','T','C','G']:
                for letter4 in ['A','T','C','G']:
                    if stemLength > 4: # Do stuff with this later.
                        for letter5 in ['A','T','C','G']: 
                            if stemLength > 5:
                                for letter6 in ['A','T','C','G']:
                                    if stemLength > 6:
                                        for letter7 in ['A','T','C','G']:    
                                            gcRatio = (int((letter1=='C' or letter1=='G'))+int((letter2=='C' or letter2=='G'))+int((letter3=='C' or letter3=='G'))+(letter4=='C' or letter4=='G')+int((letter5=='C' or letter5=='G'))+(letter6=='C' or letter6=='G')+(letter7=='C' or letter7=='G'))/7.0
                                            stem = letter1 + letter2 + letter3 + letter4 + letter5 + letter6 + letter7
                                            arrayOfStems.append(stem)
                                    else: 
                                        letter7 = ""
                                        if 0.7 <= gcRatio and 0.8 >= gcRatio:
                                            gcRatio = (int((letter1=='C' or letter1=='G'))+int((letter2=='C' or letter2=='G'))+int((letter3=='C' or letter3=='G'))+(letter4=='C' or letter4=='G')+int((letter5=='C' or letter5=='G'))+(letter6=='C' or letter6=='G'))/6.0     
                                            stem = letter1 + letter2 + letter3 + letter4 + letter5 + letter6 + letter7
                                            arrayOfStems.append(stem)
                            else:
                                letter6 = ""
                                letter7 = ""
                                gcRatio = (int((letter1=='C' or letter1=='G'))+int((letter2=='C' or letter2=='G'))+int((letter3=='C' or letter3=='G'))+(letter4=='C' or letter4=='G')+int((letter5=='C' or letter5=='G')))/5.0     
                                if 0.7 <= gcRatio and 0.8 >= gcRatio:
                                    stem = letter1 + letter2 + letter3 + letter4 + letter5 + letter6 + letter7
                                    arrayOfStems.append(stem)
                    else:
                        letter5 = ""
                        letter6 = ""
                        letter7 = ""
                        gcRatio = (int((letter1=='C' or letter1=='G'))+int((letter2=='C' or letter2=='G'))+int((letter3=='C' or letter3=='G'))+(letter4=='C' or letter4=='G')+int((letter5=='C' or letter5=='G')))/4.0     
                        if 0.7 <= gcRatio and 0.8 >= gcRatio:
                            stem = letter1 + letter2 + letter3 + letter4 + letter5 + letter6 + letter7
                            arrayOfStems.append(stem)
    return(arrayOfStems)

# Returns array of the complements for all the stems. 
def createAntiStems(arrayOfStems):
    revComplimentKey = arrayOfStems[0].maketrans("ACTG", "TGAC")
    antiStems = []
    for stem in arrayOfStems:
        antiStems.append(stem.translate(revComplimentKey)[::-1])
    return antiStems

arrayOfStems = createStems(stemLength)
arrayOfAntiStems = createAntiStems(arrayOfStems)

# String containing all possible beacon sequences.
allPossibleBeacons = ''
allBeaconsList = []

# Populate allPossibleBeacons; separate each sequence with a semicolon.
# I do this because it's the required input format of QuikFold. 
for index in range(0, len(arrayOfStems)):
    allBeaconsList.append(arrayOfStems[index] + probeSequence + arrayOfAntiStems[index])
    allPossibleBeacons = allPossibleBeacons + arrayOfStems[index] + probeSequence + arrayOfAntiStems[index] + ';'

lengthOfBeacon = len(probeSequence) + 2 * stemLength 

# Go to the Quikfold server website where we'll fold all possible beacons.     
driver.get("http://unafold.rna.albany.edu/?q=DINAMelt/Quickfold")
multiSequenceBox = driver.find_element_by_name("seq") 
multiSequenceBox.send_keys(allPossibleBeacons) # Enter sequence into box
temperatureBox = driver.find_element_by_name('temp')
temperatureBox.send_keys(Keys.BACKSPACE + Keys.BACKSPACE + str(temperature)) # Enter folding temperature 
driver.find_element_by_xpath('//*[@id="f"]/p[9]/input[1]').click() # Submit job
driver.find_element_by_xpath('//*[@id="content"]/div/div[2]/div/div/p[4]/a[1]').click() # Get ct file

# Output of the DNA folding software QuikFold
foldOutput = driver.find_element_by_xpath('/html/body/pre').text

# Class to make foldOutput nicer... eventually.
class mfoldSequence():
    # Sequence is the sequence of a given probe. Currently it's a list of chars. 
    # Hybridizations is an array of the indices of all hybridized bases in the sequence.
    # DG is a float representing the Gibbs Free Energy of the sequence.
    def __init__(self,sequence,hybridizations,DG):
        self.sequence = sequence
        self.hybridizations = hybridizations
        self.DG = DG

# Array to hold all mfold output sequence objects
mfoldSequences = []

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

def getHybridizations(foldOutput,beaconIndex):
    indicesOfHybridizations = []
    for rowIndex in range(0,len(foldOutput[beaconIndex])):
        if int(foldOutput[beaconIndex][rowIndex][4]) != 0:
            indicesOfHybridizations.append(rowIndex)
    return(indicesOfHybridizations)

foldOutput, dGs = splitString(foldOutput)

# Need to adjust this loop to only add the lowest-DG duplicate of a given sequence. 
# Or not, because hey look at that, I made a function to do it. 
for index in range(0, len(foldOutput)):
    mfoldSequences.append(mfoldSequence(list(np.transpose(foldOutput[index])[1]), getHybridizations(foldOutput, index), dGs[index]))

"""
# Remove all duplicate sequences within mfoldSequences. Keep the ones with minimum dG. Thankfully, they're already sorted by dG. 
for i in range(0, len(mfoldSequences)-1):
    #isoSequences = mfoldSequences[i] Was this line ever necessary? I don't think so. Bring it back if the whole thing breaks.
    for j in range(i+1, len(mfoldSequences)):
        # RUH ROH: The below line creates an index out of bound error. 
        if mfoldSequences[i].sequence == mfoldSequences[j].sequence and mfoldSequences[i].DG < mfoldSequences[j].DG: # THIS LINE produces an out of range error.
            mfoldSequences.remove(mfoldSequences[j])
"""

# Remove duplicate sequences; retain the minimum DG sequence if dupes are removed.
# This only works if the sequence list passed is sorted with minimum DGs first; this
# is the case in our instance (yay!). 
def deDupe(sequences):
    deDuped = []
    for i in range(0, len(sequences)-1): 
        if sequences[i].sequence not in deDuped: 
            deDuped.append(sequences[i]) 
    return(deDuped) 

mfoldSequences = deDupe(mfoldSequences)

# Takes the .ct output from Quikfold and determines which sequences are acceptable molecular beacons. 
# Does this by making sure that only the stem nucleotides are hybridized. 
# BAD????????????????????????????????????????????????????????????????????
def selectForGoodBeacons(foldOutput):
    goodBeaconIndices = []
    for beaconIndex in range(0,len(mfoldSequences)):
        indicesOfHybridizations = mfoldSequences[beaconIndex].hybridizations
        correctNumberOfHybridizations = (len(indicesOfHybridizations) == 2 * stemLength)
        correctIndicesOfHybridizations = (indicesOfHybridizations == list(range(stemLength)) + list(range(lengthOfBeacon-stemLength, lengthOfBeacon))) # um what the falcon is this
        if correctNumberOfHybridizations and correctIndicesOfHybridizations: 
            # Make sure your logic is correct about the output counting a hybridization for each member of hybridized pairs.
            goodBeaconIndices.append(beaconIndex)
    return(goodBeaconIndices)
        
# Indices of foldOutput which correspond to sequences within foldOutput with
# the required folding structure. 
goodBeaconIndices = selectForGoodBeacons(foldOutput) 

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

bestSequence = mfoldSequences[bestBeaconIndices[0]].sequence

# does this help? I am dumb and can't print a string that makes me happy
copyOfBestSequence = bestSequence
bestSequence = ""
for character in copyOfBestSequence:
    bestSequence = bestSequence + character
    
print(bestSequence) 

# Close browser and driver by user input. 
closeinator = input("Hit enter to close. \n")

# HEY YOU. Read this. 
# If there is a runtime error, run these methods in the console. 
driver.close() # Close the Chrome window and terminate the driver. 
driver.quit()  # Probably overkill but I don't like background processes. 

"""
# CHECK OUTPUT: 
for sequence in mfoldSequences:
    print(sequence.hybridizations)
    
getHybridizations(foldOutput,5)
"""