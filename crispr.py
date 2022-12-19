'''
Command line before starting:
sudo apt install python3
sudo apt install python3-is-python
sudo apt install python3-pip
pip install selenium
pip install numpy
sudo apt install ncbi-blast+
pip install networkx
pip install matplotlib
'''
import os
from os.path import exists
from datetime import datetime
import time
import numpy as np
import concurrent.futures
import math
import sys
from collections import defaultdict
import networkx as nx
'''
Also imported:
- Within scrape_crisprDB():
   from selenium import webdriver
   from selenium.webdriver.support.ui import Select
   from selenium.webdriver.common.by import By
   from selenium.webdriver.support.wait import WebDriverWait
   from selenium.webdriver.support import expected_conditions as EC
'''

start = time.perf_counter()

#        database search, parse db_file to fasta, make BLAST database, blast search
steps = [True,            True,                   True,                True]
 
filename = ""
filenames = {}
searchTerm = ""
output_filename = ""
inputVars = {"-s","-f"}
inputFiles = {"-o":"output","-out":"output","-fasta":"fasta","-blast_output":"blast"}
commandArgs = {"-db":0,"-fasta":1}

if(len(sys.argv) == 2):
    searchTerm = sys.argv[1]
else:
    for x in range(1,len(sys.argv)):
        if(sys.argv[x][0] == "-"):
            if sys.argv[x] in commandArgs:
                stepIndex = commandArgs[sys.argv[x]]
                for y in range(stepIndex,len(steps)):
                    steps[y] = False
            if sys.argv[x] in inputVars:
                if(sys.argv[x] == "-s"): searchTerm = sys.argv[x+1]
                elif(sys.argv[x] == "-f"): filename = sys.argv[x+1]
                x += 1
            elif sys.argv[x] in inputFiles:
                filenames[inputFiles[sys.argv[x]]] = sys.argv[x+1]
                x += 1

if(filename == ""): filename = searchTerm.replace(" ", "-").lower()   

if(not steps[0]):
    steps[0] = True
    from selenium import webdriver # imports selenium, which is needed for web interaction
    from selenium.webdriver.support.ui import Select
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.wait import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC
    now = str(datetime.now())
    
        #initializes chrome webdriver
    options = webdriver.ChromeOptions()
    options.add_experimental_option('excludeSwitches', ['enable-logging'])
    driver = webdriver.Chrome(options=options)
    driver.get("https://crisprcas.i2bc.paris-saclay.fr/MainDb/straindb")
    driver.implicitly_wait(30)
    
        # sets search conditions
    driver.find_element(By.ID, 'db-list-filter').send_keys(searchTerm)
    driver.find_element(By.CLASS_NAME, "btn-group").click()
    driver.find_element(By.XPATH, "//*[@id=\"renderbody\"]/div/div/div[1]/div[1]/div[3]/span/div/ul/li[6]/a/label").click()
    driver.find_element(By.XPATH, "//body").click()
    
    wait = WebDriverWait(driver, 10)
    
        #gets the number at the bottom of the page displaying the # of results
    wait.until(EC.text_to_be_present_in_element((By.ID, "strain-dt_info"), "filtered"))
    entries = driver.find_element(By.ID, "strain-dt_info").text
    entries = entries[entries.index("of")+3:]
    entries = entries[:entries.index(" ")]
    entries = int(entries)
    
    if(entries == 0):
        print("No results found.")
        quit()

        #grab each strain and get its internal reference in the database by parsing the javascript call
        #used for my initial purposes, but this needs to be reworked a little to pick something better than xpath to reference the elements by, otherwise can cause errors
    entryList = []
    for entry in range(entries):
        assembly = driver.find_element(By.XPATH, "/html/body/div/div/div/div[1]/div[2]/div[1]/div[2]/table/tbody/tr[" + str(entry+1) + "]/td[2]/div/a")
        html = assembly.get_attribute("innerHTML")
        html = html[html.index("(")+2:html.index(")")-1]
        entryList.append(html)
    
    def read_strainHTML(strainHTML):
        strainHTML = strainHTML[strainHTML.index("<h3=")+7:]
        strain = strainHTML[:strainHTML.index("</div>")].strip().replace("\n", " ")
        strain = strain.replace("  ", " ")
        
        strainHTML = strainHTML[strainHTML.index("Evidence<br>level"):]
        crisprHTMLs = []
        while("<td>CRISPR</td>" in strainHTML):
            strainHTML = strainHTML[strainHTML.index("<td>CRISPR</td>") + 4:]
            crisprHTMLs.append(strainHTML[strainHTML.index("<td>"):strainHTML.index("</tr>")])
        crisprList = [strain]
        for crispr in crisprHTMLs:
            if(int(crispr[crispr.index("<td id=\"elCol\">")+15:crispr.index("<td id=\"elCol\">")+16]) == 4):
                crispr = crispr[crispr.index("onclick=\"elSelected('")+21:]
                crisprID = crispr[crispr.index(">")+1:crispr.index("<")]
                crispr = crispr[:crispr.index("'")]
                crisprList.append([crisprID, crispr])
        return crisprList
    
    def read_crisprHTML(crisprHTML):
        crisprHTML = crisprHTML[crisprHTML.index("DR Consensus"):]
        crisprHTML = crisprHTML[crisprHTML.index("seq"):]
        consensus = crisprHTML[crisprHTML.index(">")+1:crisprHTML.index("<")].strip()
    
        spacerList = []
        while("<div>>spacer" in crisprHTML):
            crisprHTML = crisprHTML[crisprHTML.index("<div>>spacer") + 5:]
            crisprHTML = crisprHTML[crisprHTML.index("<div>")+5:]
            spacerList.append(crisprHTML[:crisprHTML.index("<")])
            crisprHTML = crisprHTML[crisprHTML.index("<div>>spacer") + 5:]
    
        return [consensus, spacerList]
    
    straindb = []
    for entry in entryList:
        driver.get("https://crisprcas.i2bc.paris-saclay.fr/MainDb/StrainData?id=" + entry)
        driver.implicitly_wait(10)
        html = driver.page_source
        strainOutput = read_strainHTML(html)
        crisprArrays = [strainOutput[0]]
        for x in range(1,len(strainOutput)):
            driver.get("https://crisprcas.i2bc.paris-saclay.fr/MainDb/ElementCrisprData?id=" + strainOutput[x][1])
            driver.implicitly_wait(10)
            html = driver.page_source
            crisprOutput = read_crisprHTML(html)
            crisprArrays.append([strainOutput[x][0],crisprOutput[0],crisprOutput[1]])
        straindb.append(crisprArrays)
    
    db_scrapeFile = searchTerm.replace(" ", "-").lower() + ".db_scrape"
    
    f = open('./Files/' + db_scrapeFile, "w")
    f.write(now)
    f.write(str(straindb))
    f.close()
    
    print("Search successful. Raw data may be found at ./Files/" + db_scrapeFile)

f = open('./Files/' + filename + ".db_scrape")
db_scrape = f.read()
f.close()

#reads the file output for the database scrape, it's just an array > string printed saved to a text file. So this code parses the string back into a useful array & throws it into the Strain object.
straindb = {}
separatorsdb = {}
crisprs = {}
crisprStrains = {}
spacerSeqs = {}
flag = False
while(bool("']]]" in db_scrape) & flag == False):
    db_scrape = db_scrape[db_scrape.index("['")+2:]
    strain = db_scrape[:db_scrape.index("'")]
    strainCrisprs = []
    while(db_scrape.index("'") < db_scrape.index("']]]")):
        db_scrape = db_scrape[db_scrape.index("['")+2:]
        crisprID = db_scrape[:db_scrape.index("'")]
        db_scrape = db_scrape[db_scrape.index("', ")+4:]
        separatorsdb[crisprID] = db_scrape[:db_scrape.index("'")]
        spacers = set()
        while(db_scrape.index("'") < db_scrape.index("']]")):
            db_scrape = db_scrape[db_scrape.index("'")+1:]
            db_scrape = db_scrape[db_scrape.index("'")+1:]
            spacerSequence = (db_scrape[:db_scrape.index("'")])
            spacerName = str(len(straindb)) + "_" + crisprID.replace("_","-") + "_" + str(len(spacers))
            spacerSeqs[spacerName] = spacerSequence
            spacers.add(spacerName)
            db_scrape = db_scrape[db_scrape.index("'"):]
        crisprs[crisprID] = spacers
        crisprStrains[crisprID] = strain
        strainCrisprs.append(crisprID)
        if("[" not in db_scrape):
            flag = True
            break
        
    straindb[strain] = np.array(strainCrisprs)

#takes the Strain objects and writes them into a fasta file for the database, with each spacer as a separate entry
if(not steps[1]): 
    filetext = ""
    for strain in straindb:
        for crispr in straindb[strain]:
            spacers = crisprs[crispr]
            for spacer in spacers:
                filetext += ">" + spacer + "\n" + spacerSeqs[spacer] + "\n"
    filetext = filetext[:len(filetext)-1]

    f = open(filenames.setdefault("fasta",'./Files/' + filename + ".fasta"), "w")
    f.write(filetext)
    f.close()


#code to execute blast searches as on the command line. 
if(not steps[3]):
    if(not steps[2]):
        os.system("makeblastdb -in " + filenames.setdefault("fasta",'Files/' + filename + ".fasta") + " -parse_seqids -title \"" + filename + "\" -dbtype nucl > Files/catch.txt")
        while not exists(filenames.setdefault("fasta",'./Files/' + filename + ".fasta") + ".ndb"):
            time.sleep(1)
        os.remove("Files/catch.txt")
    evalue = 10
    os.system("blastn -query " + filenames.setdefault("fasta",'./Files/' + filename + ".fasta") +" -db " + filenames.setdefault("fasta",'./Files/' + filename + ".fasta") + " -evalue " + str(evalue) + " -outfmt \"6 qseqid sseqid qlen slen sstart send evalue mismatch\" -out " + filenames.setdefault("blast",'./Files/' + filename + "-blastn.txt"))
    while not exists(filenames.setdefault("blast",'./Files/' + filename + "-blastn.txt")):
        time.sleep(1)
    f = open(filenames.setdefault("blast",'./Files/' + filename + "-blastn.txt"))
    blastOutput = str(f.read())
    f.close()
    f = open(filenames.setdefault("blast",'./Files/' + filename + "-blastn.txt"),"w")
    f.write("qseqid\tsseqid\tqlen\tslen\tsstart\tsend\tevalue\tmismatch\n" + blastOutput)
    f.close()
else:
    f = open(filenames.setdefault("blast", 'Files/' + filename +"-blastn.txt"))
    blastOutput = str(f.read())
    f.close()
    blastOutput = blastOutput[blastOutput.index("mismatch")+9:]

#blast output to list
lines = blastOutput.split("\n")
blastOutput = []
for x in range(len(lines)):
    intermediate = lines[x].split("\t")
    if(len(intermediate) == 8):
        if(not intermediate[0] == intermediate[1]):
            blastOutput.append(intermediate)
del lines
del intermediate


matchlist = []
seen = set()
for line in blastOutput:
    if(line[0] in seen):
        pass
    else:
        if(line[1] in seen):
            for setN in range(len(matchlist)):
                if(line[1] in matchlist[setN]):
                    matchlist[setN].add(line[0])
                    seen.add(line[0])
                    break                        
        else:
            matchlist.append({line[0], line[1]})
            seen.update({line[0], line[1]})
blastMap = matchlist
""" condensedList = []
while(len(blastMap) > 0):
    connections = blastMap.pop(0)
    searched = set()
    while(len(searched) < len(connections)):
        spacer = next(iter(connections-searched))
        index = len(blastMap) -1
        while(index > 0):
            if(any(spacer == newSpacer for newSpacer in blastMap[index])):
                connections.update(blastMap[index])
                blastMap.remove(blastMap[index])
            index += -1
        searched.add(spacer)
    condensedList.append(connections)
blastMap = condensedList """

print(len(blastMap))

spreadsheetStrains = defaultdict(set)
spreadsheetCrisprs = defaultdict(set)
spreadsheetSpacers = {}
for line in blastMap:
    matches = len(line)
    for spacer in line:
        crisprID = spacer[spacer.index("_")+1:]
        crisprID = crisprID[:crisprID.index("_")].replace("-","_")
        spreadsheetStrains[crisprStrains[crisprID]].add(crisprID)
        spreadsheetCrisprs[crisprID].add(spacer)
        spreadsheetSpacers[spacer] = np.array([spacerSeqs[spacer],matches])

filetext = "Strain Name\tCRISPR Array\tSpacer Sequence\tNumber of Matches\tSpacers with Matches\tTotal Spacers"

for strain in sorted(spreadsheetStrains):
    filetext += "\n" + strain + "\t\t\t\t"
    matchCount = 0
    totalCount = 0
    for crispr in spreadsheetStrains[strain]:
        matchCount += len(spreadsheetCrisprs[crispr])
        totalCount += len(crisprs[crispr])
    filetext += str(matchCount) + "\t" + str(totalCount)
    for crispr in sorted(spreadsheetStrains[strain]):
        filetext += "\n\t" + crispr + "\t\t\t" + str(len(spreadsheetCrisprs[crispr])) + "\t" + str(len(crisprs[crispr]))
        for spacer in sorted(spreadsheetCrisprs[crispr]):
            filetext += "\n\t\t" + spreadsheetSpacers[spacer][0] + "\t" + str(spreadsheetSpacers[spacer][1])

f = open(filenames.setdefault("output",'./Files/' + filename + "-out.txt"), "w")
f.write(filetext)
f.close()

def strain(spacerID):
    intermediate = spacerID[spacerID.index("_")+1:]
    intermediate = intermediate[:intermediate.index("_")].replace("-","_")
    return crisprStrains[intermediate]


pairs = {}
for line in matchlist:
    while(len(line) > 1):
        strain1 = strain(line.pop())
        for spacer in line:
            pair = frozenset([strain1, strain(spacer)])
            if(len(pair) == 2): 
                pairs[pair] = pairs.setdefault(pair, 0)+1

#for pair in pairs: print(pair)

maxMatches = max(pairs.values())

edges = []
for pair in pairs:
    edge = tuple(pair) + ({"weight":2*pairs[pair]/maxMatches},)
    edges.append(edge)

G = nx.Graph()
G.add_nodes_from(straindb.keys())
G.add_edges_from(edges)

nx.write_graphml(G, "Files/" + filename + ".graphml")




finish = time.perf_counter()
print(finish-start)
print("done!")
