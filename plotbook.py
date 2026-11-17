# plotbook.py
#
#   Plots to .pdf based on input file.
#
#   Syntax:
#   python plotbook.py [num threads] [input file] [file with list of out files to process]
#
#   example:
#   python ..\..\..\scripts\plotbook.py 20 plotbook.txt list.txt
#
#   ----------------------------------------------------------------------------
#   Schema for the plotbook input file (comma-separated):
#
#   1:  Type, see dictPlotTypes for available types.
#   2:  Plot title.
#   3:  List of arguments (usually channel names).
#
#   ---- Example plotbook.txt file below:
#   PEPL,Area PE-PL (NB NS),NB,NS
#   ANGL_ALL,NS machine angles (all),MAC1
#   PELEC_PLOAD,Area PElec (NB NS),NB,NS
#   ANGL_M,Thermal mach angl,TC3,TR6,PA1,LG3
#   MACH_MW,Thermal mach (P),TC3,TR6,PA1,LG3
#   MACH_MVAR,Thermal mach (Q),TC3,TR6,PA1,LG3
#   FREQ,Frequency,103H_345,103H_138,103H_69,104H_138,67N_345,67N_230,91H_138,101S_345,410N_345
#   BUSVOLT,Voltage,103H_345,103H_138,103H_69,104H_138,67N_345,67N_230,91H_138,101S_345,410N_345
#   MACH_MW,Wind IRs (P),574,597
#   MACH_MVAR,Wind IRs (Q),574,597
#   BUSVOLT,Wind IRs (V),574,597
#   ---- End of example
#
#   Revision history: ----------------------------------------------------------
#
#   2026/09/28  PVH
#   Modified for py 3.11
#
#   2022/12/21  PVH
#   Moved to using custom tools library.
#
#   2020/06/18  PVH
#   Added dynamic marker offsets (still needs work).
#
#   2020/06/16  PVH
#   Added IR542 to wind machine filter list.
#
#   2020/06/15  PVH
#   Added filtering option for ANGL_ALL plots.
#
#   2020/06/08  PVH
#   Initial release.
#

# ---- Start of imports -------------------------------------------------------
#  System start here
import os, sys, socket
import csv          # csv file processing
import re           # regex
import glob         # pathname pattern matching
import pdb          # python debugger
import logging      # python logging
import math         # math functions

import datetime, time         # timekeeping

# multiproc
from concurrent.futures import ProcessPoolExecutor

# More plotting imports
import numpy as np      # Array manipulations
import pandas as pd     # Data analysis library
import random           # Random number generator
import gc               # Garbage collection
import itertools        # Cycling library
from matplotlib.backends.backend_pdf import PdfPages    # .pdf support
import matplotlib.ticker as ticker                      # Axis tick control
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---- Initialize program settings for PSSE run outside of GUI ----------------
# Hostnames
nameWstn        = "EMA1902003L"                             # Workstation name

# Install paths
locPsseDef      = "C:\Program Files (x86)\PTI\PSSE34\PSSPY27"
locPsseWksn     = "C:\Program Files (x86)\pti\psse34\PSSPY27"
nameHost = socket.gethostname()

# Add  PSSEBIN dir, decision is made on where to source additional tools
if nameHost == nameWstn:
    PSSE_LOCATION = locPsseWksn
    # Include Peter"s tools path
    sys.path.append("c:/work/resources/planning/psse/psse_api")
# If not on Peter"s computer, assume they"re packed locally
else:
    PSSE_LOCATION = locPsseDef

sys.path.append(PSSE_LOCATION)
os.environ["PATH"] = os.environ["PATH"] + ";" +  PSSE_LOCATION

# ---- Resuming imports -------------------------------------------------------
import psse35
import psspy
import dyntools
import redirect

redirect.psse2py()
psspy.throwPsseExceptions = True

# Include custom tools
import ierrRef
import custTools

# ---- End of   imports -------------------------------------------------------


# ---- Start of global definitions --------------------------------------------
outPdfDir       = "outpdf"                      # output pdf directory

dictPlotTypes   = { "ANGL_ALL"  : "ANGL_ALL",
                    "ANGL_M"    : "ANGL_M",
                    "CORRIDOR"  : "CORRIDOR",
                    "LINE_MW"   : "LINE_MW",
                    "LINE_MVAR" : "LINE_MVAR",
                    "BUSVOLT"   : "BUSVOLT",
                    "FREQ"      : "FREQ",
                    "PEPL"      : "PEPL",
                    "PMECH"     : "PMECH",
                    "PELEC"     : "PELEC",
                    "PLOAD"     : "PLOAD",
                    "PELEC_PLOAD" : "PELEC_PLOAD",
                    "CMLD_ALL"  : "CMLD_ALL",
                    "MACH_MW"   : "MACH_MW",
                    "MACH_MVAR" : "MACH_MVAR",

                    # "Hidden" ones
                    "NBNS_TIE"  : "NBNS_TIE",
                    "NS_CBX"    : "NS_CBX",
                    "NS_ONI"    : "NS_ONI",
                    "NS_ONS"    : "NS_ONS",
                    "120H_SVC_MVAR":    "120H_SVC_MVAR",
                    
                    # Not developed yet
                    "CMLD_A"    : "CMLD_A",
                    "CMLD_B"    : "CMLD_B",
                    "CMLD_C"    : "CMLD_C",
                    "CMLD_D"    : "CMLD_D",
                    "CMLD_PE"   : "CMLD_PE",
                    "CMLD_S"    : "CMLD_S"
                    }
# ---- End of   global definitions --------------------------------------------


# ---- Start of Chart() -------------------------------------------------------
# Class for each channel output file plot, containing many subplots.
class Chart:
    def __init__(   self, itemY,
                    itemX       = None,
                    titleFig    = None,
                    titlePage   = None,
                    labelX      = None,
                    labelY      = None,
                    limX        = None,
                    limY        = None,
                    stepX       = None,
                    stepY       = None,
                    bLegend     = None,
                    type        = None
                    ):
        self.itemY      = itemY
        self.itemX      = itemX
        self.titleFig   = titleFig
        self.titlePage  = titlePage
        self.labelX     = labelX
        self.labelY     = labelY
        self.limX       = limX
        self.limY       = limY
        self.stepX      = stepX
        self.stepY      = stepY
        self.bLegend    = bLegend
        self.type       = type
# ---- End of   Chart() -------------------------------------------------------


# ---- Start of ChannelPlot() -------------------------------------------------
# Class for each channel output file, containing the data and the chart data.
class ChannelPlot:
    def __init__(   self,
                    datCh,
                    myPlot
                    ):
        self.datCh      = datCh
        self.myPlot     = myPlot                                                # Array of Chart class
# ---- End of   ChannelPlot() -------------------------------------------------


# ---- Start of retMathRegEx() ------------------------------------------------
# Returns regex substring for CMLD channel pattern.
#
#   Arguments:
#       type:   String from dictPlotTypes[]
def retMathRegEx(type):
    if(type == dictPlotTypes["CMLD_ALL"]):
        substr = "^199.* P$"
    elif(type == dictPlotTypes["CMLD_A"]):
        substr = "^199.*A P$"
    elif(type == dictPlotTypes["CMLD_B"]):
        substr = "^199.*B P$"
    elif(type == dictPlotTypes["CMLD_C"]):
        substr = "^199.*C P$"
    elif(type == dictPlotTypes["CMLD_D"]):
        substr = "^199.*D P$"
    elif(type == dictPlotTypes["CMLD_PE"]):
        substr = "^199.*E P$"
    elif(type == dictPlotTypes["CMLD_S"]):
        substr = "^199.*S P$"
    elif(type == dictPlotTypes["NBNS_TIE"]):
        substr = "^L(8001:410N|6535:410N|6536:410N).*MW$"
    elif(type == dictPlotTypes["NS_CBX"]):
        substr = "^L(8004:101S|7003:3C|7004:3C|7005:3C|6515:2C).*MW$"
    elif(type == dictPlotTypes["NS_ONI"]):
        substr = "^L(8003:67N|7003:67N|7005:67N|7019:67N|6503:1N).*MW$"
    else:
        substr = ".*"

    return substr
# ---- End of   retMathRegEx() ------------------------------------------------


# ---- Start of calcCMLDComps() -----------------------------------------------
# Returns CMLD component summations.
#
#   Arguments:
#       type:   String from dictPlotTypes[]
#       df:     Pandas dataframe with all data.
def calcCMLDComps(type, df):
    # Get regex substring
    substrCMLD      = retMathRegEx(singleplt.type)
    
    # Filter columns to match regex pattern
    #filteredCMLDdat = df.filter(regex = substrCMLD)
    
    # Add the newly calculated column to the dataframe
    #df[type]   = filteredCMLDdat.sum(axis = 1)
    
    print(df.head())
    
    return df
# ---- End of   calcCMLDComps() -----------------------------------------------


# ---- Start of readRunFile() -------------------------------------------------
# Reads plotbook file.
#   Arguments:
#       fNameRunFile:   Filename for the plot input.
#   Returns:
#       Array, with the contents of each row.
def readRunFile(fNameRunFile):
    fRunList    = open(fNameRunList, "r")
    contents    = []
    parsedRow   = []

    with open(fNameRunList) as csvFile:
        hCSV = csv.reader(csvFile, delimiter=",")

        for row in hCSV:
            parsedRow = []

            for ctr, arg in enumerate(row, start = 1):
                parsedRow.append(arg)

            contents.append(parsedRow)

    return contents
# ---- End of   readRunFile() -------------------------------------------------


# ---- Start of readOutFileList() ---------------------------------------------
# Reads plotbook file.
#   Arguments:
#       fNameRunFile:   Filename for the plot input.
#   Returns:
#       Array, with the contents of each row.
def readOutFileList(fNameOutFile):
    fOutList    = open(fNameOutFile, "r")
    contents    = []

    with open(fNameOutFile, "r") as fOutList:
        contents = [line.rstrip() for line in fOutList]

    return contents
# ---- End of   readOutFileList() ---------------------------------------------


# ---- Start of getRelevantChans() --------------------------------------------
# Reads plotbook file.
#   Arguments:
#       context:    Dictionary from dictPlotTypes.
#       datChID:    Dictionary of chan #'s & names to search.
#       details:    Array of detailed channel identifers to search with.
#   Returns:
#       List of relevant channel numbers.
def getRelevantChans(context, datChID, details):
    newDict = {}

    
    if(context == dictPlotTypes["ANGL_ALL"]):
        # Create dictionary, with only machine angles
        substr  = re.compile(".*ANGLM.*")
        tempDict = {}

        for key, value in datChID.items():
            if(substr.match(value)):
                tempDict[key] = value

        # Build machine filter (wind machines and any special exclusions)
        #strFilter = "19C|106W|109S|110W|92N|89N|91N|93N|102V|120C|98V|IR542|301NS_MLP1|301NS_MLP2|IR760|IR761|IR673|IR686|IR670|112W|IR668|IR739|IR618_742|IR728|IR732|IR826|IR677|50N_ST1|BESS|IR669|IR675"
        strFilter = "301NS_MLP"
        
        # Add on any extra exclusions to the filter pattern, if supplied
        if(len(details) > 0):
            strXtraExcl = "|".join(details)
            strFilter   = strFilter + "|" + strXtraExcl

        # Filter out wind machines, reverse lookup
        substr  = re.compile("^(?!.*(" + strFilter + ").*)")
        for key, value in tempDict.items():
            if(substr.match(value)):
                #print("Match: " + str(key) + ", " + value)
                newDict[key] = value
            #else:
            #    print("No match: " + str(key) + ", " + value)

        #print(newDict)

    elif(   (context == dictPlotTypes["ANGL_M"]) or
            (context == dictPlotTypes["CORRIDOR"]) or
            (context == dictPlotTypes["BUSVOLT"]) or
            (context == dictPlotTypes["LINE_MW"]) or
            (context == dictPlotTypes["LINE_MVAR"]) or
            (context == dictPlotTypes["FREQ"]) or
            (context == dictPlotTypes["PEPL"]) or
            (context == dictPlotTypes["PMECH"]) or
            (context == dictPlotTypes["PELEC"]) or
            (context == dictPlotTypes["PLOAD"]) or
            (context == dictPlotTypes["PELEC_PLOAD"]) or
            (context == dictPlotTypes["CMLD_ALL"]) or
            (context == dictPlotTypes["CMLD_A"]) or
            (context == dictPlotTypes["CMLD_B"]) or
            (context == dictPlotTypes["CMLD_C"]) or
            (context == dictPlotTypes["CMLD_D"]) or
            (context == dictPlotTypes["CMLD_PE"]) or
            (context == dictPlotTypes["CMLD_S"]) or
            (context == dictPlotTypes["MACH_MW"]) or
            (context == dictPlotTypes["MACH_MVAR"])
            ):
        # Create dictionary, with only class names
        if(context == dictPlotTypes["ANGL_M"]):
            substr  = re.compile("^ANGLM.*$")

        if(context == dictPlotTypes["CORRIDOR"]):
            substr  = re.compile("^.*MW$")
            
        elif(context == dictPlotTypes["BUSVOLT"]):
            substr  = re.compile("^VOLTB.*$")
            
        elif(context == dictPlotTypes["LINE_MW"]):
            substr  = re.compile("^L.*MW$")
            
        elif(context == dictPlotTypes["LINE_MVAR"]):
            substr  = re.compile("^L.*MVAR$")
            
        elif(context == dictPlotTypes["LINE_MVAR"]):
            substr  = re.compile("^L.*MVAR$")
            
        elif(context == dictPlotTypes["FREQ"]):
            substr  = re.compile("^FREQ.*$")
            
        elif(context == dictPlotTypes["PEPL"]):
            substr  = re.compile("^PE-PL.*$")
            
        elif(context == dictPlotTypes["PMECH"]):
            substr  = re.compile("^PMECH.*$")
            
        elif(context == dictPlotTypes["PELEC"]):
            substr  = re.compile("^PELEC.*$")
            
        elif(context == dictPlotTypes["PLOAD"]):
            substr  = re.compile("^PLOAD.*$")
            
        elif(context == dictPlotTypes["PELEC_PLOAD"]):
            substr  = re.compile("^(PELEC|PLOAD).*$")
            
        elif(context == dictPlotTypes["CMLD_ALL"]):
            substr  = re.compile("^.*(A P|B P|C P|D P|S P|E P)$")
            
        elif(context == dictPlotTypes["CMLD_A"]):
            substr  = re.compile("^.*A P$")
            
        elif(context == dictPlotTypes["CMLD_B"]):
            substr  = re.compile("^.*B P$")
            
        elif(context == dictPlotTypes["CMLD_C"]):
            substr  = re.compile("^.*C P$")
            
        elif(context == dictPlotTypes["CMLD_D"]):
            substr  = re.compile("^.*D P$")
            
        elif(context == dictPlotTypes["CMLD_PE"]):
            substr  = re.compile("^.*E P$")
            
        elif(context == dictPlotTypes["CMLD_S"]):
            substr  = re.compile("^.*S P$")
            
        elif(context == dictPlotTypes["MACH_MW"]):
            substr  = re.compile("^PELEM.*$")
            
        elif(context == dictPlotTypes["MACH_MVAR"]):
            substr  = re.compile("^(QELEM|120H_SVC).*$")                        # Kludge to include the SVC

        tempDict = {}

        for key, value in datChID.items():
            if(substr.match(value)):
                tempDict[key] = value

        # Filter through the class for each reference to the specified bus names.
        for bus in details:
            substr  = re.compile(".*" + bus + ".*")

            for key, value in tempDict.items():
                if(substr.match(value)):
                    newDict[key] = value

    return newDict
# ---- End of   getRelevantChans() --------------------------------------------


# ---- Start of transformDataframe() ------------------------------------------
# Applies math to dataframe, if necessary.
def transformDataframe(type, values):

    valuesProc = np.array(values)
    
    if(type == dictPlotTypes["FREQ"]):
        valuesProc = (1 + valuesProc) * 60

    #print(valuesProc.head())
    
    return valuesProc
    
# ---- End of   transformDataframe() ------------------------------------------


# ---- Start of makePlot() ----------------------------------------------------
# Writes the plots to .pdf.
#
#   Arguments:
#   myCases:    An array of class ChannelPlot.
#   fNamePdf:   Filename, without extension, for the output pdf.
#   appLog:     Logging handler.
#   layout:     Tuple (x, y), of # subplots on each page. Experimental, haven't tested other than (2, 2).
def makePlot(myCases, fNamePdf = "myPlots", appLog = 0, layout = (2, 2)):

    # Definitions
    SIZE_SMALL      = 6                                                         # Font size, small
    SIZE_MEDIUM     = 9                                                         # Font size, medium
    SIZE_BIGGER     = 11                                                        # Font size, bigger
        
    pageSizeX       = 11                                                        # Page width, in inches
    pageSizeY       = 8.5                                                       # Page height, in inches
    paddingX        = 0.1                                                       # Subplot distance X
    paddingY        = 0.1                                                       # Subplot distance Y
    paddingXLabel   = 0                                                         # Subplot x-axis label padding
    paddingYLabel   = 2                                                         # Subplot y-axis label padding
    subplotLeft     = 0.05                                                      # Subplot positioning left
    subplotRight    = 0.99                                                      # Subplot positioning right
    subplotBottom   = 0.04                                                      # Subplot positioning bottom
    subplotTop      = 0.95                                                      # Subplot positioning top

    tickLen         = 12                                                        # Subplot axis tick length
    
    lineWidth       = 0.9                                                       # Plot line weight

    # Markers
    markerSize      = 5                                                         # Marker size
    arrMarker       = [ ".", "+", "x",   2, "*",   3,   4, "d", "|",   5,   6,  # Array of markers to cycle through
                          7, "1", "2", "o", "v", "^", "<", ">", "3", "4", "8",
                        "s", "p", "P", "h", "H", "X", "D", "_",   0,   1,   8,
                          9,  10,  11]

    # Update rc
    params = {  "font.size"         : SIZE_SMALL,                               # Controls default text sizes
                "axes.titlesize"    : SIZE_SMALL,                               # Fontsize of the axes title
                "axes.labelsize"    : SIZE_SMALL,                               # Fontsize of the x and y labels
                "axes.titlepad"     : paddingYLabel,                            # Subplot y padding
                "xtick.labelsize"   : SIZE_SMALL,                               # Fontsize of the tick labels
                "ytick.labelsize"   : SIZE_SMALL,                               # Fontsize of the tick labels
                "legend.fontsize"   : SIZE_SMALL,                               # Legend fontsize
                "figure.titlesize"  : SIZE_BIGGER                               # Fontsize of the figure title
                }
    plt.rcParams.update(params)

    with PdfPages(fNamePdf + ".pdf") as hPdf:
        fig         = plt.figure()
        ctrPg       = 0                                                         # Global page counter
        ctrPgCase   = 0                                                         # Page counter for each case
        ctrPgInst   = 0                                                         # Subplot instance counter (per page)

        pdfCaseLim  = 8                                                         # Max cases to pdf file
        totPkgCase  = len(myCases)

        # pdf metadata
        hMeta   = hPdf.infodict()
        hMeta["Title"]          = "Title"
        hMeta["Author"]         = "Author"
        hMeta["Subject"]        = "Subject"
        hMeta["CreationDate"]   = datetime.datetime.today()
        hMeta["ModDate"]        = datetime.datetime.today()

        # Figure setup
        fig.set_size_inches([pageSizeX, pageSizeY])

        # Roll through each case
        for ctrCase, myPlots in enumerate(myCases):
            ctrPgCase       = 0
            ctrPdf          = 0
            
            totCasePages    = int(math.ceil(float(len(myPlots.myPlot)) / float(layout[0] * layout[1]))) # Amount of pages per case, counting fr 0
            totPdfPages     = (totCasePages) * len(myCases)                                 # Total amount of pages in plotbook, counting fr 0
            
            #print("totCasePages: " + str(totCasePages) + ", len(myPlots.myPlot): " + str(len(myPlots.myPlot)) + ", layouts: " + str(layout[0] * layout[1]) + ".")
            #appLog.info("Plotting case " + str(ctrCase + 1) + "/" + str(totPkgCase) + ".")

            # Roll through the plots in each case
            for ctr, singleplt in enumerate(myPlots.myPlot):
                # Start each page off right, with a some general information
                if(ctrPgInst == 0):
                    # Title, upper-left
                    fig.text(   x           = 0.02,
                                y           = 0.98,
                                s           = singleplt.titlePage + " (" + str(ctrPgCase + 1) + "/" + str(totCasePages) + ")",
                                ha          = 'left',
                                fontsize    = SIZE_MEDIUM
                                )
                    # Page number/total pages, lower-right
                    fig.text(   x           = 0.99,
                                y           = 0.01,
                                s           = str(ctrPg + 1) + "/" + str(totPdfPages),
                                ha          = 'right',
                                fontsize    = SIZE_SMALL
                                )
                                
                    ctrPg       = ctrPg + 1

                plt.subplots_adjust(    left    = subplotLeft,
                                        right   = subplotRight,
                                        bottom  = subplotBottom,
                                        top     = subplotTop,
                                        wspace  = paddingX,
                                        hspace  = paddingY
                                        )

                # Create subplot
                fig.add_subplot(layout[0],
                                layout[1],
                                ctrPgInst + 1,
                                )

                # Build dataframe, starting with the time column as the index
                #df = pd.DataFrame({ "time"      : myPlots.datCh["time"]
                #                    })

                if(0):
                    df = pd.DataFrame.from_dict({"time"  : myPlots.datCh["time"]})
                    df = df.set_index("time")

                    # Add channels to the dataframe
                    for key, value in singleplt.itemY.items():
                        df[value] = myPlots.datCh[key]
                else:
                    # Build dataframe in a single operation to avoid fragmentation

                    chanData = {
                        value: myPlots.datCh[key]
                        for key, value in singleplt.itemY.items()
                    }

                    df = pd.DataFrame(
                        chanData,
                        index=myPlots.datCh["time"]
                    )

                    df.index.name = "time"

                # Set current axis
                ax  = plt.gca()                            
                ax.set_title(singleplt.titleFig)

                # Additional customizations, like tick spacing
                # Bus voltage plots
                if(singleplt.type == dictPlotTypes["BUSVOLT"]):
                    plt.ylim(0.0, 1.2)                                          # Fix y-axis range
                    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))        # Set minor tick spacing
                
                # Frequency plots
                elif(singleplt.type == dictPlotTypes["FREQ"]):
                    plt.ylim(58, 62)                                            # Fix y-axis range
                    ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))        # Set major tick spacing
                    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))       # Set minor tick spacing
                
                # MW plots (large)
                elif(singleplt.type == dictPlotTypes["CORRIDOR"]):
                    ax.yaxis.set_minor_locator(plt.MultipleLocator(100))        # Set minor tick spacing

                # MW plots (line flow)
                elif(singleplt.type == dictPlotTypes["LINE_MW"]):
                    ax.yaxis.set_minor_locator(plt.MultipleLocator(20))        # Set minor tick spacing

                # MW plots (CMLD)
                elif(singleplt.type == dictPlotTypes["CMLD_ALL"]):
                    ax.yaxis.set_minor_locator(plt.MultipleLocator(50))        # Set minor tick spacing

                # Machine MW PU output plots (small)
                elif(   (singleplt.type == dictPlotTypes["MACH_MW"])
                        ):
                    plt.ylim(0, 2.5)                                            # Fix y-axis range
                    ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))        # Set major tick spacing
                    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))       # Set minor tick spacing

                # Machine MVAR PU output plots (small)
                elif(   (singleplt.type == dictPlotTypes["MACH_MVAR"])
                        ):
                    plt.ylim(-2, 3)                                             # Fix y-axis range
                    ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))        # Set major tick spacing
                    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))        # Set minor tick spacing

                # All other plots use auto ticks/y-axis ranges
                else:
                    ax.relim()                                                  # Get axis limits
                    ax.autoscale()                                              # Reset scales
                    
                #plt.xlim(   xmin    = 0.5,
                #            xmax    = 3.5
                #            )
                    
                ax.xaxis.set_minor_locator(plt.MultipleLocator(1))              # Set x-axis minor tick spacing

                # Columns math, if necessary
                # CMLD summations
                if(singleplt.type == dictPlotTypes["CMLD_ALL"]):

                    # For each CMLD component:
                    #   Get regex pattern.
                    #   Filter the main dataframe to match the regex pattern.
                    #   Perform a summation of all channels.
                    substrCMLD      = retMathRegEx(dictPlotTypes["CMLD_ALL"])
                    allCMLDChans    = df.filter(regex = substrCMLD)
                    df[dictPlotTypes["CMLD_ALL"]]   = allCMLDChans.sum(axis = 1)
                    df[dictPlotTypes["CMLD_ALL"]]   = df[dictPlotTypes["CMLD_ALL"]].apply(lambda x: x * (100))   # Convert to MW
                    
                    substrCMLD      = retMathRegEx(dictPlotTypes["CMLD_A"])
                    allCMLDA        = df.filter(regex = substrCMLD)
                    df[dictPlotTypes["CMLD_A"]]     = allCMLDA.sum(axis = 1)
                    df[dictPlotTypes["CMLD_A"]]   = df[dictPlotTypes["CMLD_A"]].apply(lambda x: x * (100))   # Convert to MW
                    
                    substrCMLD      = retMathRegEx(dictPlotTypes["CMLD_B"])
                    allCMLDB        = df.filter(regex = substrCMLD)
                    df[dictPlotTypes["CMLD_B"]]     = allCMLDB.sum(axis = 1)
                    df[dictPlotTypes["CMLD_B"]]   = df[dictPlotTypes["CMLD_B"]].apply(lambda x: x * (100))   # Convert to MW
                    
                    substrCMLD      = retMathRegEx(dictPlotTypes["CMLD_C"])
                    allCMLDC        = df.filter(regex = substrCMLD)
                    df[dictPlotTypes["CMLD_C"]]     = allCMLDC.sum(axis = 1)
                    df[dictPlotTypes["CMLD_C"]]   = df[dictPlotTypes["CMLD_C"]].apply(lambda x: x * (100))   # Convert to MW
                    
                    substrCMLD      = retMathRegEx(dictPlotTypes["CMLD_D"])
                    allCMLDD        = df.filter(regex = substrCMLD)
                    df[dictPlotTypes["CMLD_D"]]     = allCMLDD.sum(axis = 1)
                    df[dictPlotTypes["CMLD_D"]]   = df[dictPlotTypes["CMLD_D"]].apply(lambda x: x * (100))   # Convert to MW
                    
                    substrCMLD      = retMathRegEx(dictPlotTypes["CMLD_PE"])
                    allCMLDPE       = df.filter(regex = substrCMLD)
                    df[dictPlotTypes["CMLD_PE"]]     = allCMLDPE.sum(axis = 1)
                    df[dictPlotTypes["CMLD_PE"]]   = df[dictPlotTypes["CMLD_PE"]].apply(lambda x: x * (100))   # Convert to MW
                    
                    substrCMLD      = retMathRegEx(dictPlotTypes["CMLD_S"])
                    allCMLDS        = df.filter(regex = substrCMLD)
                    df[dictPlotTypes["CMLD_S"]]     = allCMLDS.sum(axis = 1)
                    df[dictPlotTypes["CMLD_S"]]   = df[dictPlotTypes["CMLD_S"]].apply(lambda x: x * (100))   # Convert to MW

                    # Remove all other columns
                    df = df.loc[:, df.columns.intersection([dictPlotTypes["CMLD_ALL"],
                                                            dictPlotTypes["CMLD_A"],
                                                            dictPlotTypes["CMLD_B"],
                                                            dictPlotTypes["CMLD_C"],
                                                            dictPlotTypes["CMLD_D"],
                                                            dictPlotTypes["CMLD_PE"],
                                                            dictPlotTypes["CMLD_S"]
                                                            ])]

                # Corridor summations for:
                #   NB/NS tie
                #   CBX
                #   ONI
                elif(singleplt.type == dictPlotTypes["CORRIDOR"]):
                    # For each corridor component:
                    #   Get regex pattern.
                    #   Filter the main dataframe to match the regex pattern.
                    #   Perform a summation of all channels.
                    substrCorridor  = retMathRegEx(dictPlotTypes["NBNS_TIE"])
                    allChansTie     = df.filter(regex = substrCorridor)
                    df[dictPlotTypes["NBNS_TIE"]]   = allChansTie.sum(axis = 1)
                    df[dictPlotTypes["NBNS_TIE"]]   = df[dictPlotTypes["NBNS_TIE"]].apply(lambda x: x * (-1))   # Invert, so tie value is from NS perspective
                    
                    substrCorridor  = retMathRegEx(dictPlotTypes["NS_CBX"])
                    allChansCBX     = df.filter(regex = substrCorridor)
                    df[dictPlotTypes["NS_CBX"]]   = allChansCBX.sum(axis = 1)
                    
                    substrCorridor  = retMathRegEx(dictPlotTypes["NS_ONI"])
                    allChansONI     = df.filter(regex = substrCorridor)
                    df[dictPlotTypes["NS_ONI"]]   = allChansONI.sum(axis = 1)
                    df[dictPlotTypes["NS_ONI"]]   = df[dictPlotTypes["NS_ONI"]].apply(lambda x: x * (-1))   # Invert, so ONI is using same convention as SCADA
                    
                    #print(df.head())

                    # Remove all other columns
                    df = df.loc[:, df.columns.intersection([dictPlotTypes["NBNS_TIE"],
                                                            dictPlotTypes["NS_CBX"],
                                                            dictPlotTypes["NS_ONI"]
                                                            ])]

                # Convert SVC MVAR output from system base to PU:
#                elif(singleplt.type == dictPlotTypes["MACH_MVAR"]):
#                    df[dictPlotTypes["120H_SVC_MVAR"]]   = df[dictPlotTypes["120H_SVC_MVAR"]].apply(lambda x: x / 100)   # Convert from SBASE 100 to PU

                # End of columns math

                if(1):
                    ctrMarker = 0

                    # Experimental: Plot all series by iterating through the dataframe
                    for (colName, colData) in df.items():
                        #if(colName != "time"):

                        # Show marker or not, depending on legend visibility
                        if(singleplt.bLegend == True):
                            markerUse = arrMarker[ctrMarker]
                            ctrMarker = ctrMarker + 1
                        else:
                            markerUse = "None"

                        # Apply math, if necessary
                        if(singleplt.type == dictPlotTypes["FREQ"]):
                            df[colName] = (1 + df[colName]) * 60

                        # Calculate x-axis marker offset
                        index   = df.index
                        numRows = len(index)
                        
                        offsetSpacing   = int(0.058 * numRows)
                        offsetOffset    = int(0.029 * numRows)
                        
                        markEvery       = offsetSpacing + (ctrMarker * offsetOffset)
                        
                        #print("numRows:       " + str(numRows))
                        #print("offsetSpacing: " + str(offsetSpacing))
                        #print("offsetOffset:  " + str(offsetOffset))
                        #print("mark every:    " + str(markEvery))
                        #print("\n")

                        df.plot(    kind        = "line",
                                    y           = colName,
                                    ax          = ax,
                                    linewidth   = lineWidth,
                                    legend      = singleplt.bLegend,
                                    marker      = markerUse,
                                    markersize  = markerSize,
                                    markevery   = markEvery     # To add some offset to markers
                                    )

                        # Reset if ctrMarker exceeds marker array size
                        if(ctrMarker >= len(arrMarker)):
                            ctrMarker = 0
                else:
                    # Just plot the dataframe
                    df.plot(    kind        = "line",
                                x           = "time",
                                ax          = ax,
                                linewidth   = 1,
                                legend      = singleplt.bLegend
                                )

                # Plot formatting
                plt.xlabel(singleplt.labelX, labelpad = paddingXLabel)
                plt.ylabel(singleplt.labelY, labelpad = paddingYLabel)
                
                ax.tick_params( axis        = "x",
                                which       = "both",
                                direction   = "in",
                                top         = True,
                                bottom      = True,
                                length      = tickLen
                                )
                ax.tick_params( axis        = "x",
                                which       = "minor",
                                length      = tickLen / 2
                                )
                ax.tick_params( axis        = "y",
                                which       = "both",
                                direction   = "in",
                                left        = True,
                                right       = True,
                                length      = tickLen
                                )
                ax.tick_params( axis        = "y",
                                which       = "minor",
                                length      = tickLen / 2
                                )

                # Set x-axis range
                #ax.set_xlim([0.5, 3.5])

                # Re-adjust legend to order the series' are drawn
                if(singleplt.bLegend == True):
                   # Get current handles and labels
                    current_handles, current_labels = plt.gca().get_legend_handles_labels()

                    # Sort or reorder the labels and handles
                    reversed_handles    = reversed(current_handles)
                    reversed_labels     = reversed(current_labels)

                    # call plt.legend() with the new values
                    #plt.legend(reversed_handles, reversed_labels)
                    plt.legend(current_handles, current_labels)

                # Check if we've filled up the page and create a new if we have.
                if((ctrPgInst >= ((layout[0] * layout[1]) - 1))):
                    ctrPgInst = 0

                    # Only add a new page if there's more subplots to be added.
                    if(ctr <= (len(myPlots.myPlot) - 2)):
                        # Add figure to .pdf
                        hPdf.savefig(fig)

                        # Clear figure for next figure
                        plt.clf()

                    # Update page counters
                    ctrPgCase   = ctrPgCase + 1
                else:
                    ctrPgInst = ctrPgInst + 1


            #plt.show(fig)

            # Add figure to .pdf
            hPdf.savefig(fig)
            plt.clf()
            
            # Reset page counters
            ctrPgInst   = 0
            ctrPdf      = ctrPdf + 1

        # Cleanup
        #hPdf.close()
# ---- End of   makePlot() ----------------------------------------------------


# ---- Start of processOutFile() ----------------------------------------------
def processOutFile(args):

    (
        ctrOut,
        fNameOut,
        totCases,
        argPlot
    ) = args

    print(f"\n----- out file: {fNameOut}")

    myPlots = []

    objOut = dyntools.CHNF(fNameOut, outvrsn=0)

    datShortTitle, datChID, datCh = objOut.get_data(
        outfile=fNameOut
    )

    for arg in argPlot:

        showLegend = True

        if(arg[0] == dictPlotTypes["ANGL_ALL"]):
            listUse = arg[2:]
            chanList = getRelevantChans(
                arg[0],
                datChID,
                listUse
            )

            labelY = "Angle (deg)"
            showLegend = False

        elif(
            (arg[0] == dictPlotTypes["ANGL_M"]) or
            (arg[0] == dictPlotTypes["BUSVOLT"]) or
            (arg[0] == dictPlotTypes["LINE_MW"]) or
            (arg[0] == dictPlotTypes["LINE_MVAR"]) or
            (arg[0] == dictPlotTypes["FREQ"]) or
            (arg[0] == dictPlotTypes["PEPL"]) or
            (arg[0] == dictPlotTypes["PMECH"]) or
            (arg[0] == dictPlotTypes["PELEC"]) or
            (arg[0] == dictPlotTypes["PLOAD"]) or
            (arg[0] == dictPlotTypes["PELEC_PLOAD"]) or
            (arg[0] == dictPlotTypes["CMLD_ALL"]) or
            (arg[0] == dictPlotTypes["CORRIDOR"]) or
            (arg[0] == dictPlotTypes["MACH_MW"]) or
            (arg[0] == dictPlotTypes["MACH_MVAR"])
        ):

            listUse = arg[2:]
            chanList = getRelevantChans(
                arg[0],
                datChID,
                listUse
            )

            labelMap = {
                "ANGL_M": "Angle (deg)",
                "BUSVOLT": "V PU",
                "LINE_MW": "MW",
                "LINE_MVAR": "MVAR",
                "FREQ": "Freq (Hz)",
                "PEPL": "PE-PL (MW)",
                "PMECH": "PMech",
                "PELEC": "MW",
                "PLOAD": "MW",
                "PELEC_PLOAD": "MW",
                "CMLD_ALL": "MW",
                "CORRIDOR": "MW",
                "MACH_MW": "MW (PU)",
                "MACH_MVAR": "MVAR (PU)"
            }

            labelY = labelMap[arg[0]]

        sTitleUse = arg[1]
        sTitlePageUse = fNameOut[:-4]

        myPlots.append(
            Chart(
                titleFig=sTitleUse,
                titlePage=sTitlePageUse,
                itemY=chanList,
                labelX="Time (s)",
                labelY=labelY,
                bLegend=showLegend,
                type=arg[0]
            )
        )

    myCase = ChannelPlot(
        datCh=datCh,
        myPlot=myPlots
    )

    outBaseName = os.path.splitext(
        os.path.basename(fNameOut)
    )[0]

    pdfName = os.path.join(
        outPdfDir,
        f"plt{ctrOut:04d}_{outBaseName}"
    )

    makePlot(
        [myCase],
        fNamePdf=pdfName
    )

    return pdfName
# ---- End of   processOutFile() ----------------------------------------------


# ---- Start of main() --------------------------------------------------------
# def main:
if __name__ == "__main__":

    bContentsFound  = 0                                                         # Flag for contents found.
    ctrPdfs         = 1                                                         # Running counter of output pdfs
    fNamePdf        = "myPlots"                                                 # Output pdf filename, w/o extension

    # create output directory
    if not os.path.exists(outPdfDir):
        os.makedirs(outPdfDir)

    # Benchmarking time placeholders
    timeStart   = time.time()
    timeEnd     = 0
    
    # Init session log
    formatterLog = logging.Formatter('%(asctime)s %(levelname)s %(funcName)s(%(lineno)d) %(message)s', datefmt='%y-%m-%d %H:%M:%S')

    fNameSession    = "plotbook_" + datetime.datetime.now().strftime("%y%m%d_%H%M%S")     #plotbook_yymmdd_hhmmss.txt

    bLog2File       = 0
    # Set up file handler
    if(bLog2File):
        hLog = logging.FileHandler(fNameSession + custTools.dictExt["ext_txt"])
        hLog.setFormatter(formatterLog)
        hLog.setLevel(logging.INFO)

    # Set up console stream handler
    sLog = logging.StreamHandler()
    sLog.setFormatter(formatterLog)
    sLog.setLevel(logging.INFO)

    # Get logger
    appLog = logging.getLogger('root')
    appLog.setLevel(logging.INFO)

    # Add both handlers
    if(bLog2File):
        appLog.addHandler(hLog)
    appLog.addHandler(sLog)

    # Parse arguments
    # Get number of threads
    if(len(sys.argv) >= 3):
        numThreads      = int(sys.argv[1])

    # Process plot file
    if(len(sys.argv) >= 3):
        fNameRunList    = sys.argv[2]
        bContentsFound  = 1
        
    # Process file list, if specified
    if(len(sys.argv) >= 4):
        fNameOutList    = sys.argv[3]
        bFileListFound  = 1
        
        fNamePdf        = sys.argv[3].rsplit(".", 1)[0]             # Remove extension (after last ".")

    # Process file list (either default everything or specific list
    if(bFileListFound):
        listFilesOut    = readOutFileList(fNameOutList)

    # Or just process all .out files in the directory
    else:
        listFilesOut    = glob.glob("*.out")

    if(bContentsFound == 1):
        argPlot = readRunFile(fNameRunList)

        totCases = len(listFilesOut)

        # For each channel file
        jobArgs = []

        for ctrOut, fNameOut in enumerate(
                listFilesOut,
                start=1):

            jobArgs.append(
                (   ctrOut,
                    fNameOut,
                    totCases,
                    argPlot
                )
            )

        appLog.info(f"Starting {numThreads} worker processes.")

        with ProcessPoolExecutor(
                max_workers=numThreads) as pool:

            for pdfName in pool.map(
                    processOutFile,
                    jobArgs,
                    chunksize=5):

                appLog.info(
                    f"Completed {pdfName}.pdf"
                )

    # Calculate end time and duration
    timeEnd     = time.time()
    timeDiff    = (timeEnd - timeStart)
    
    appLog.info("Plots complete. Duration was " + str("{:.0f}".format(timeDiff)) + " seconds.")
# ---- End of main() ----------------------------------------------------------
