#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code author: Katie E. Blise
Original Date: March, 2021
Revised Date: September, 2021

This .py file contains all code needed to reproduce the figures in the paper:
"Single-Cell Spatial Proteomics Analyses of Head and Neck Squamous Cell Carcinoma Reveal Tumor Heterogeneity and Immune Architectures Associated with Clinical Outcome"


This .py file contains multiple functions:
    getCsvs() - generates a list of all 47 tumor regions analyzed, in correct order for figure generation
    cellClassRename - generates dictionary with gated cell labels converted to biologically relevant cell phenotypes
    densityToCsv() - create csvs with counts of cell types for all tumor regions, to generate figure 1
    calculateMixing() - calculates mixing score for all tumor regions, to generate figure 3
    calculateCold() - calculates which tumor regions are labeled 'cold,' to generage figure 3
    pdCounts() - create csvs with counts of PD-1+ and PD-L1+ cells for all tumor regions, to generate figure 3
    kiCounts() - create csvs with counts of Ki-67+ cells for all tumor regions, to generate figure 3
    makeNeighborhoods() - create csv with cellular neighborhood compositions for each seed cell within distance threshold, to generate figure 4
    elbowMethod() - gerate elbow plot to determine optimal number of clusters for k-means clustering
    clusterNeighborhoods() - create csv with results from k-means clustering on cellular neighborhoods, to generate figure 4
    clusterCountPerROI() - create csvs with counts of cellular neighborhood clusters for all tumor regions, to generate figure 4
    fig1() - generates Figures 1b-g
    fig2() - generates Figures 2a-d
    fig3() - generates Figures 3b-i
    fig4() - generates Figures 4b-d, f-g
    table3() - generates Table 3
    supp1() - generates Supplemental Figure 1c
    supp2() - generates Supplemental Figures 2a-e
    supp3() - generates Supplemental Figures 3b-f
    

Note: This program assumes the following items live in the same directory as this .py file:
    -'dfCreated' folder
    -'figures' folder
    -'data' folder, which houses: all mIHC data (.csv files) and clinical data file (.csv)  

This program is intended for Python version 3.
"""



def getCsvs():
    
    '''
    This function generates a list containing all 47 tumor regions analyzed.
    Note: all tumor regions are assumed to be downloaded and stored in the 'data' folder.
    
    Input parameters:
        None
    
    Output:
        Returns csvList = a list that contains all 47 tumor regions analyzed, in order of pt#, pre/post, roi
    '''
    
    #import packages
    import os
    
    #get list of all files in current directory/data
    cwd = os.getcwd()
    files = os.listdir(cwd+'/data')
    
    #keep only tumor regions in csvList
    filesUpdate = []

    for file in files:
        if 'roi' in file:
            filesUpdate.append(file[0:-4]) #add file, drop '.csv' for figure code
        
    #create csvList to store correct csvs, in correct order for further analyses
    csvList = []
    
    treatList = ['pre','post']
    
    #iterate through 9 patients
    for i in range(1,10):
        
        pt = 'pt'+str(i)
        
        for treat in treatList:
            
            #create pt# string to match on
            region = pt+treat
            
            toAddList = [roi for roi in filesUpdate if region in roi]
            toAddList.sort()
            
            #add csvs in order to csvList
            csvList = csvList + toAddList
    
    return csvList



def cellClassRename():
    
    '''
    This function returns a dicionary of matched cell class IDs (A/B/C/etc) with their corresponding cell class names (CD8+ T cell/CD4+ T-helper cell/B cell/etc).
    
    Input parameters:
        None        
        
    Output:
        returns nameDict = a dictionary of {cell class IDs:cell class names} 
    '''
    
    nameDict = {'A':'CD8+ T Cell','B':'CD4+ T Helper','C':'B Cell','D':'Macrophage','E':'Other Immune','F':'Other Non-Immune','G':'Noise','H':'Neoplastic Tumor','J':'Granulocyte','K':'CD4+ Regulatory T Cell','N':'⍺SMA+ Mesenchymal','X':'Antigen Presenting Cell'}
       
    return nameDict



def densityToCsv(cellList):

    '''
    This function calculates the density of cells present in a given list of csvs. It saves the results to csvs - both for all regions and averaged across regions (for a single tumor).

    Input parameters:
        cellList = list of cell classes to count; eg ['A','B','C','H']

    Outputs:
        Saves two csvs to the dfCreated folder; one with density of all regions and one with density averaged across regions
    '''

    #import packages
    import pandas as pd
    import os

    #get csvList
    csvList = getCsvs()

    #get current directory
    cwd = os.getcwd()

    ##create new df and add clinical data to it
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv', index_col=0)
    dfClin.set_index('sample',inplace=True) #set index labels to be the patient column

    #create a new df to store clinical data and counts
    dfDensity = pd.DataFrame()
    dfDensity['patient'] = csvList
    dfDensity.set_index('patient',inplace=True) #set index labels to be the patient column
    dfDensity['dtr'] = dfClin['dtr']
    dfDensity['dtr'] = pd.to_numeric(dfDensity['dtr'],errors='coerce')

    #set up to store densities in value's list of countDict
    countDict = {}
    for cell in cellList:
        countDict[cell] = [] #empty list to hold cell counts for each cell type we care about

    #loop through csvs and calculate counts
    for file in csvList:
        df = pd.read_csv(cwd+'/data/'+file+'.csv',index_col=0)

        #get area of the region
        area = dfClin.loc[file,'area']

        #only include specified cell types in the df
        df = df[[cell in cellList for cell in df['class']]]

        #get counts of all cells in the df
        freq = df['class'].value_counts()

        for cell in cellList: #only include cells that we want

            if cell in freq: #check if the cell type exists in the file
                countDict[cell].append(freq[cell]/area) #add the frequency/area (density) of the cell type to the dictionary value's list
            else:
                countDict[cell].append(0) #otherwise add zero as the count

    #put counts into dfDensity
    for c,countList in countDict.items():
        dfDensity[c] = countList

    #save to csv for all regions
    dfDensity.to_csv(cwd+'/dfCreated/dfDensity_all.csv')

    #create new dataframe with averaged numbers for every ROI within one patient
    dfDensity.index = dfDensity.index.str.slice(0,6) #update patient column with just the first six characters in string; so you can use groupby

    #groupby same patient (diff ROI) and take the average of the corresponding rows
    dfDensityAvg = dfDensity.groupby(['patient']).mean()

    #save dfAvg to csv
    dfDensityAvg.to_csv(cwd+'/dfCreated/dfDensity_avg.csv')



def calculateMixing(distThresh,cellList):
    
    '''
    This function calculates the mixing score for all tumor regions and saves the results to a csv in the 'dfCreated' folder.
    It assumes that interactions are counted between cells defined in the cellList parameter and tumor cells (H).
    Note that the raw counts of cell interactions is doubled, but when doing the mixing score calculation, the division negates this.       

    Input parameters:
        distThresh = distance threshold to count cells as interacting; radius; in pixels (2px = 1µm)
        cellList = list of cell types to include as cells counted in interaction tallies in the mixing score calculation; for all immune: ['A','B','C','D','E','J','K','X']
    
    Outputs:
        Saves two csvs to the dfCreated folder; one with mixing scores of all regions and one with scores averaged across regions  
    '''
    
    #import packages
    import pandas as pd
    import os
    import numpy as np
    from scipy import spatial

    #get csvList
    csvList = getCsvs()
    
    #get current directory
    cwd = os.getcwd()
        
    ##create new df and add clinical data to it
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv', index_col=0)
    
    #empty lists to store clinical data from dfClin
    dtrList = []
    
    #create a new df to store clinical data and scores
    dfMix = pd.DataFrame()
    dfMix['patient'] = csvList
    dfMix.set_index('patient',inplace=True) #set index labels to be the patient column
    
    #loop through each csv in csvList to get clinical values of only the csvs included in the list
    for pt in csvList:
        idx = dfClin.loc[dfClin['sample'] == pt].index[0]
        dtrList.append(dfClin['dtr'].iloc[idx]) #add to dtr list
    
    #put lists into dfScore dataframe and convert all values to int and any missing values to NaN values
    dfMix['dtr'] = dtrList
    dfMix['dtr'] = pd.to_numeric(dfMix['dtr'],errors='coerce')    
    
    #calculate mix scores - store values to list to add to dfMix
    mixList = []
    
    #loop through each tumor region in csvList and calculate its mixing score
    for file in csvList:

        print('Calculating mixing score for '+file+'...')
        
        #read original tumor region csv
        df = pd.read_csv(cwd+'/data/'+file+'.csv', index_col=0)
        
        #create list of cell types to keep in the df; will include the specified cells from cellList parameter plus the tumor population
        keepCells = cellList + ['H']
            
        #create filtered dataframe using keepCells list
        filt_df = df[df['class'].isin(keepCells)]
        
        #create np array of just x,y coordinates
        ptsArray = filt_df[['Location_Center_X','Location_Center_Y']].values
        
        #create kdtree
        tree = spatial.KDTree(ptsArray)
        
        #create empty dictionary to hold index:list of neighboring cells' indeces
        neighDict = {}
        
        #determine cells within radius using KDtree query_ball_point
        for i in range(len(ptsArray)):
            neighs = tree.query_ball_point(ptsArray[i],distThresh) #get neighbors within distance
            if i in neighs:
                neighs.remove(i) #don't include itself as a neighbor
            neighDict[i] = neighs #index:list of neighboring cells' indeces
        
    
        #tumor cells are type H; cellList is specified as a parameter    
        tumList = ['H'] #all runs have H as the tumor cell type
            
        #set counts to zero to start
        tumImmCount = 0
        immImmCount = 0
 
        #calculate mixing score
        for key,val in neighDict.items(): #key = index of cell1; val = list of neighboring cells' indeces
            
            class1 = filt_df['class'].values[key] #get original classification of cell1
    
            for n in val: #loop through each neighbor indeces
                class2 = filt_df['class'].values[n] #get original classification of the neighboring cells
    
                if (class1 in tumList and class2 in cellList) or (class2 in tumList and class1 in cellList): #tumor:immune interaction
                    tumImmCount = tumImmCount + 1
        
                elif class1 in cellList and class2 in cellList: #immune:immune interaction
                    immImmCount = immImmCount + 1
    
        #avoid dividing by zero
        if immImmCount != 0:
            mixScore = tumImmCount/immImmCount
    
        elif immImmCount == 0:
            mixScore = np.nan #assign NaN value to avoid dividing by zero; so it's properly evaluated (ignored) in the mean calculation
                    
        #add mixing score to mixList
        mixList.append(mixScore)

    #put mixList into dfMix
    dfMix['mixScore'] = mixList
    
    #save dfMix to csv; so when you restart kernel you don't have to rerun analysis on all tumors
    dfMix.to_csv(cwd+'/dfCreated/dfMixingScores_all.csv')
    
    #create new dataframe with averaged numbers for every ROI within one patient
    dfMix.index = dfMix.index.str.slice(0,6) #update patient column with just the first six characters in string; so you can use groupby
    
    #groupby same patient (diff ROI) and take the average of the corresponding rows
    dfMixAvg = dfMix.groupby(['patient']).mean()
    
    #save dfMixAvg to csv
    dfMixAvg.to_csv(cwd+'/dfCreated/dfMixingScores_avg.csv')



def calculateCold():

    '''
    This function calculates which tumor regions are designated cold and saves new csvs that exclude the cold regions.
    It uses the same threshold (ratio) as that in Keren et al.
    
    Input parameters:
        None
    
    Outputs:
        Saves two csvs to the dfCreated folder; one with mixing scores of all regions and one with scores averaged across regions, although both exclude cold regions
    '''
    
    #import packages
    import pandas as pd
    import os
    
    #get csvList
    csvList = getCsvs()
    
    #get current directory
    cwd = os.getcwd()
    
    #set immune cell types to count up
    immList = ['A','B','C','D','E','J','K','X']
    
    immCountDict = {} #to store file: immune cell density
    
    #get area values to calculate immune density
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv', index_col=0)
    dfClin.set_index('sample',inplace=True) #set index labels to be the sample column
    
    for file in csvList:
    
        #get original roi file
        df = pd.read_csv(cwd+'/data/'+file+'.csv', index_col=0)
    
        #get corresponding area of that roi
        area = dfClin.loc[file,'area']
        
        #get raw counts of each class in the roi
        counts = df['class'].value_counts()
        
        #reset for each csv; tally total immune cells
        count = 0
        
        for i in immList:
            if i in counts:
                count = count + counts[i]    
        #add immune count/area (density) to dictionary for each roi
        immCountDict[file] = count/area   
    
    #identify csvs with a count of fewer than 390.625 immune cells; these are the 'cold' regions; Keren et al used <250 immune per 0.64mm2 area which is equal to <390.625
    excludeList = [] #empty list to store cold regions
    for key,value in immCountDict.items():
        if value < 390.625: #threshold determined using Keren et al's threshold (adjusting for image size)
            #print(key)
            excludeList.append(key) #add region to the exclude list
    
    print('Regions to exclude:',excludeList)        
    
    #call mixing score df
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv',index_col=0)
    
    #get rid of cold ROIs
    dfFilt = dfMix.loc[~dfMix.index.isin(excludeList)]
    
    ##save cold excluded results to new csvs
    dfFilt.to_csv(cwd+'/dfCreated/dfMixingScoresExcludeCold_all.csv')
    
    #reset index to then use groupby - for average
    dfFilt.index = dfFilt.index.str.slice(0,6)
    
    #groupby same patient and average
    dfFiltAvg = dfFilt.groupby(['patient']).mean()
    
    #save updated dfFiltAvg averaged across patients to a new csv
    dfFiltAvg.to_csv(cwd+'/dfCreated/dfMixingScoresExcludeCold_avg.csv')



def pdCounts():
    
    '''
    This function counts the number of cells in each region that are gated positive for PD-1 and PD-L1.
    It calculates percentage counts: for example (PD-1+ CD8 T cells)/(all CD8 T cells)
    Results for % counts are saved to csvs - one for all regions and one averaged across regions (per tumor).
    
    Input parameters:
        None
        
    Outputs:
        Saves two csvs to the dfCreated folder; one with % PD-1+ and PD-L1+ counts of all regions and one with % counts averaged across regions    
    '''
    
    import pandas as pd
    import os
    
    #get csvList
    csvList = getCsvs()
    
    #get current directory
    cwd = os.getcwd()
     
    #get clinical df
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv', index_col=0)
    
    #empty lists to store clinical data from dfClin
    dtrList = []
    
    #create a new df to store clinical data and % counts
    dfCount = pd.DataFrame()
    dfCount['patient'] = csvList
    dfCount.set_index('patient',inplace=True) #set index labels to be the patient column
    
    #loop through each csv in csvList to get clinical values of only the csvs included in the list
    for pt in csvList:
        idx = dfClin.loc[dfClin['sample'] == pt].index[0]
        dtrList.append(dfClin['dtr'].iloc[idx]) #add to dtr list
    
    #put lists into dfScore dataframe and convert all values to int and any missing values to NaN values
    dfCount['dtr'] = dtrList
    dfCount['dtr'] = pd.to_numeric(dfCount['dtr'],errors='coerce')
    
    #set up to store pd1/pdl1 raw counts in value's list of countDict
    countDict = {}
    dfRaw = pd.DataFrame() #empty df to store raw counts
    dfRaw['patient'] = csvList
    dfRaw.set_index('patient',inplace=True) #set index labels to be the patient column
    
    #list all cell types to count pd-1 and pd-l1 expression
    cellList = ['A','B','C','D','E','J','K','X','N','H']   
    
    #create empty lists to hold cell counts for each cell type/PD-status we care about
    for cell in cellList:
        countDict[cell+'_pd1+'] = [] #PD-1+
        countDict[cell+'_pdl1+'] = [] #PD-L1+
    
    #loop through csvs and calculate counts
    for file in csvList:    
        #get original csv
        df = pd.read_csv(cwd+'/data/'+file+'.csv',index_col=0)
    
        for cell in cellList:
            #make smaller df of just the cell type we want
            dfCell = df[(df['class'] == cell)] #note that if cell type does not exist in the df, zeros will automatically be added as the counts
    
            #only count the cells positive for one (or both) of these markers
            #get counts of PD-1
            pd1Pos = len(dfCell[(dfCell['Cellsp_PD1p'] == 1)])
    
            #get counts of PD-L1
            pdl1Pos = len(dfCell[(dfCell['Cellsp_PDL1p'] == 1)])
        
            #NOTE: if a cell is positive for BOTH PD-1 and PD-L1 it will get counted for BOTH markers.
            #add counts to countDict
            countDict[cell+'_pd1+'].append(pd1Pos) #PD-1+
            countDict[cell+'_pdl1+'].append(pdl1Pos) #PD-L1+
    
    #put counts into dfRaw
    for c,countList in countDict.items():
        dfRaw[c] = countList    
    
    #get list of regions
    ptList = list(dfRaw.index)
    
    #get list of columns
    colList = list(dfRaw.columns)
    
    #set up to store % counts in value's list of countDict
    countDict = {}
    for c in colList:
        countDict[c] = [] #empty list to hold cell counts for each cell type/PD-status we care about
            
    #get original region's csv to get cell counts from
    for pt in ptList:
        #print(pt) #to keep track of progress
        df = pd.read_csv(cwd+'/data/'+pt+'.csv',index_col=0)
    
        #get array of counts for each cell type in the original ROI
        cellCounts = df['class'].value_counts()
            
        #calculate % of each cell
        for c in colList: #loop through each column of the dfRaw file
            
            if c[0] in cellCounts:
                countDict[c].append((dfRaw.loc[pt,c])/(cellCounts[c[0]])) #append the % count = raw count of cell type/total of that cell type
    
            else:
                countDict[c].append(0) #if cell type is not present in the original ROI, then add zero percent to the countDict for that cell type    
        
    #put counts into dfCount
    for c,countList in countDict.items():
        dfCount[c] = countList
    
    #save dfCount to csv
    dfCount.to_csv(cwd+'/dfCreated/dfPD1PDL1CountsPerc_all.csv')
    
    #create new dataframe with averaged numbers for every region within one patient
    dfCount.index = dfCount.index.str.slice(0,6)
    
    #groupby same patient (diff region) and take the average of the corresponding rows
    dfCountAvg = dfCount.groupby(['patient']).mean()
    
    #save dfCountAvg to csv
    dfCountAvg.to_csv(cwd+'/dfCreated/dfPD1PDL1CountsPerc_avg.csv')
    


def kiCounts():

    '''
    This function counts the number of cells in each region that are gated positive for Ki-67.
    It calculates percentage counts: for example (Ki-67+ CD8 T cells)/(all CD8 T cells)
    Results for % counts are saved to csvs - one for all regions and one averaged across regions (per tumor).

    Input parameters:
        None

    Outputs:
        Saves two csvs to the dfCreated folder; one with % Ki-67+ counts of all regions and one with % counts averaged across regions    
    '''

    import pandas as pd
    import os

    #get csvList
    csvList = getCsvs()

    #get current directory
    cwd = os.getcwd()

    #get clinical df data
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv', index_col=0)

    #empty lists to store clinical data from dfClin
    dtrList = []

    #create a new df to store clinical data and % counts
    dfCount = pd.DataFrame()
    dfCount['patient'] = csvList
    dfCount.set_index('patient',inplace=True) #set index labels to be the patient column

    #loop through each csv in csvList to get clinical values of only the csvs included in the list
    for pt in csvList:
        idx = dfClin.loc[dfClin['sample'] == pt].index[0]
        dtrList.append(dfClin['dtr'].iloc[idx]) #add to dtr list

    #put lists into dfScore dataframe and convert all values to int and any missing values to NaN values
    dfCount['dtr'] = dtrList
    dfCount['dtr'] = pd.to_numeric(dfCount['dtr'],errors='coerce')

    #set up to store ki67 raw counts in value's list of countDict
    countDict = {}
    dfRaw = pd.DataFrame() #empty df to store raw counts
    dfRaw['patient'] = csvList
    dfRaw.set_index('patient',inplace=True) #set index labels to be the patient column

    #list all cell types to count ki-67 expression
    cellList = ['A','B','C','D','E','J','K','X','N','H']   

    #create empty lists to hold cell counts for each cell type
    for cell in cellList:
        countDict[cell+'_ki67+'] = []

    #loop through csvs and calculate counts
    for file in csvList:
        #get original csv
        df = pd.read_csv(cwd+'/data/'+file+'.csv',index_col=0)

        for cell in cellList:
            #make smaller df of just the cell type we want
            dfCell = df[(df['class'] == cell)] #note that if cell type does not exist in the df, zeros will automatically be added as the counts

            #only count the cells positive for ki67
            ki67Pos = len(dfCell[(dfCell['Cellsp_KI67p'] == 1)])

            #add counts to countDict
            countDict[cell+'_ki67+'].append(ki67Pos)

    #put counts into dfRaw
    for c,countList in countDict.items():
        dfRaw[c] = countList    

    #get list of regions
    ptList = list(dfRaw.index)

    #get list of columns
    colList = list(dfRaw.columns)

    #set up to store % counts in value's list of countDict
    countDict = {}
    for c in colList:
        countDict[c] = [] #empty list to hold cell counts for each cell type

    #get original region's csv to get cell counts from
    for pt in ptList:
        df = pd.read_csv(cwd+'/data/'+pt+'.csv',index_col=0)

        #get array of counts for each cell type in the original ROI
        cellCounts = df['class'].value_counts()

        #calculate % of each cell
        for c in colList: #loop through each column of the dfRaw file

            if c[0] in cellCounts:
                countDict[c].append((dfRaw.loc[pt,c])/(cellCounts[c[0]])) #raw count of cell type/total of that cell type

            else:
                countDict[c].append(0) #if cell type is not present in the original ROI, then add zero percent to the countDict for that cell type    

    #put counts into dfCount
    for c,countList in countDict.items():
        dfCount[c] = countList

    #save dfCount to csv
    dfCount.to_csv(cwd+'/dfCreated/dfKI67CountsPerc_all.csv')

    #create new dataframe with averaged numbers for every region within one patient
    dfCount.index = dfCount.index.str.slice(0,6)

    #groupby same patient (diff region) and take the average of the corresponding rows
    dfCountAvg = dfCount.groupby(['patient']).mean()

    #save dfCountAvg to csv
    dfCountAvg.to_csv(cwd+'/dfCreated/dfKI67CountsPerc_avg.csv')



def makeNeighborhoods(seed,distThresh):

    '''
    This function counts the neighboring cell types within X distance of the seed cells and saves these counts to a csv.
    Counts are saved both as a raw count and as a percentage.    

    Input parameters:
        seed = seed cell type to create neighborhoods for (ex. H)
        distThresh = radius (in pixels; 2 px = 1µm) to draw around seed cell to identify neighbors in

    Output:
        One csv is saved in the 'dfCreated/' folder within the original path.
        It is saved with the following naming convention: dfNeighborhoodCluster'+seed+str(distThresh)+'.csv'
        This csv is what should be used to cluster on.
        Each row is one seed cell and the columns contain the number of neighbors it has of each cell type.
    '''

    import pandas as pd
    from scipy import spatial
    import os

    #get csvList
    csvList = getCsvs()
    
    #get current directory
    cwd = os.getcwd()
    
    #empty list to hold dictionaries to create new rows of dfClust; outside of for file in csvList loop
    allNeighList = []

    #loop through each ROI in csvList
    for file in csvList:

        #read df according to path
        df = pd.read_csv(cwd+'/data/'+file+'.csv', index_col=0)

        #create filtered dataframe without noise or 'other cells'
        filt_df = df[((df['class'] != 'G') & (df['class'] != 'F'))] #G=noise, F=other non-immune cells

        #get all possible class values; for later counting
        classOptions = sorted(list(filt_df['class'].unique()))

        #generate count variable names for final created df's columns
        classCountNames = []
        for c in classOptions: #loop through all possible neighboring cells
            classCountNames.append('count'+c)

        ##get nearest neighbors of seed cells defined by seed param
        #create np array of just x,y coordinates
        ptsArray = filt_df[['Location_Center_X','Location_Center_Y']].values

        #create kdtree
        tree = spatial.KDTree(ptsArray)
        #print('tree created.')

        #loop through each cell and check its neighbors if it's the right seed
        for i in range(len(ptsArray)):

            classType = filt_df['class'].values[i]
            #print('\n',i,classType)

            #only check neighbors if the seed cell is the desired classType
            if classType == seed:
                neighs = tree.query_ball_point(ptsArray[i], distThresh) #get neighbors within distance

                #print(neighs)
                if i in neighs:
                    neighs.remove(i) #don't include itself as a neighbor

                neighClassList = [] #empty list to store neighboring cells' classes
                 
                #loop through each neighbor and get its class; add it to neighClassList
                for j in neighs:
                    #get its class
                    neighClass = filt_df['class'].values[j]
                    neighClassList.append(neighClass)

                #get counts of neighboring cell types
                classCounts = []
                for c in classOptions: #loop through all possible neighboring cells
                    count = neighClassList.count(c) #count the specified cell type
                    classCounts.append(count) #add the counts of the specified cell type to classCounts

                #reset dictionary to hold counts of neighbors; one per seed
                seedDict = {}

                #add counts to a dictionary (one per seed cell); also add original seed cell's idx value (filt_df.iloc's index)
                seedDict['file'] = file
                seedDict['index'] = i

                #for each possible neighboring cell class, add its count to seedDict both raw and as a %
                for n in range(len(classCountNames)):
                    seedDict[classCountNames[n]] = classCounts[n] #raw count

                    if sum(classCounts) != 0:
                        seedDict[classCountNames[n]+'%'] = classCounts[n]/sum(classCounts) #percentage
                    else: #avoid division by zero if there are no neighbors
                        seedDict[classCountNames[n]+'%'] = 0 #set % to zero if there are no neighbors

                #add each seed's neighbor dictionary to the overall list; one dictionary per row of df
                allNeighList.append(seedDict)

    #create one new df to hold data for clustering; pass in allNeighList as the data; format is one row per seed cell across all csvs
    #column names from a seedDict's keys (all seedDicts have the same keys)
    dfClust = pd.DataFrame(data = allNeighList, columns = list(seedDict.keys()))
    #NOTE: only using the columns (cells present) from the final csv that was analyzed. If a cell is not present at all in this csv then its column will not be present in the final csv created.

    #convert any NaN values to zeros [note that NaN values arise when a csv lacks any of a cell type that does exist in other csvs]
    dfClust = dfClust.fillna(0)

    #store dfClust as a csv
    dfClust.to_csv(cwd+'/dfCreated/dfNeighborhoodCluster'+seed+str(distThresh)+'.csv')
    print('Csv created. Check columns to ensure all expected cells are present as columns.')



def elbowMethod(file,steps):
    
    '''
    This function runs the Elbow Method to determine the optimal number of clusters for k-means clustering.
    Input parameters:
        file = name of file to run clustering on (likely of the format dfNeighborhoodClusterH50)
        steps = max number of clusters (k) to test
     
    Output:
        saves plot to 'figures' folder
    '''
    
    import pandas as pd
    from matplotlib import pyplot as plt
    from sklearn.cluster import MiniBatchKMeans #minibatchkmeans is better when n > 10,000 samples
    import os
    
    cwd = os.getcwd()
    
    df = pd.read_csv(cwd+'/dfCreated/'+file+'.csv', index_col=0)

    #drop all rows that have no cells in the neighborhood (aka when the sum of count columns is zero)
    df['sum'] = df.iloc[:,2:].sum(axis=1)
    df = df[df['sum'] != 0]
        
    #generate column list to cluster on based on if there is a % in the column name
    colList = list(df.columns[['%' in col for col in list(df.columns)]])
    
    #get only features we want to cluster on
    df = df[colList]
    data = df.values
    
    #empty list to store error value
    wcss = []
    
    #calculate error for each k value (k=number of clusters)
    for k in range(1, steps):
        #print('running for k=',k)
    
        #generate kmeans model
        kmeans = MiniBatchKMeans(n_clusters=k, init='k-means++', max_iter=300, n_init=10, random_state=0)
    
        #fit model to data
        kmeans.fit(data)
    
        #add the sum of squares to wcss list; for plotting elbow
        wcss.append(kmeans.inertia_)
    
    #generate elbow plot
    plt.plot(range(1, steps), wcss)
    plt.title('Elbow Method')
    plt.xlabel('Number of clusters')
    plt.ylabel('WCSS')
    plt.savefig(cwd+'/figures/fig4ElbowPlot.png',format='png') #saves figure to 'figures' folder as a png file
    plt.close()
    print("Elbow plot saved to 'figures' folder.")
 


def clusterNeighborhoods(file,k):
    
    '''
    This function runs k-means clustering on a given neighborhood clustering csv.
    The results are saved to a new csv.

    Input parameters:
        file = name of file to run clustering on excluding the .csv (ex. dfNeighborhoodClusterN60)
        k = number of clusters; use elbow method to determine optimal number

    Outputs:
        One csv is saved to the 'dfCreated/' folder in the form: 'dfNeighClustered'+seedDist+'+roiType+'k'+k+'.csv'

    '''

    import pandas as pd
    from sklearn.cluster import MiniBatchKMeans
    import os

    #get current directory
    cwd = os.getcwd()

    #read csv with neighborhood data    
    df = pd.read_csv(cwd+'/dfCreated/'+file+'.csv', index_col=0)
  
    #drop all rows that have no cells in the neighborhood (aka when the sum of count columns is zero)
    df['sum'] = df.iloc[:,2:].sum(axis=1)
    dfFilt =df[df['sum'] != 0]

    #get lists from original df to add back later after clustering
    roiList = list(dfFilt['file']) #get list of ROIs in order to add to dfNoNa later
    idxList = list(dfFilt['index']) #get list of cell indices in order to add to dfNoNa later

    #generate column list to cluster on based on if there is a % in the column name; this excludes the ki67, pdl1 columns too
    colList = list(dfFilt.columns[['%' in col for col in list(dfFilt.columns)]])

    #get only features we want to cluster on
    dfFilt = dfFilt[colList]
    data = dfFilt.values

    #=k-means clustering of cells with k clusters
    kmeans = MiniBatchKMeans(n_clusters=k, init='k-means++', max_iter=300, n_init=10, random_state=0)
    predict = kmeans.fit_predict(data) #fit model to data and predict index (cluster labels); same results as first fitting and then predicting
    #print('prediction complete.')

    #add predicted cluster labels to df as a new column and see how many cells are in each cluster
    dfFilt['cluster'] = predict
    dfFilt['cluster'].value_counts() #unecessary unless you want to see how many cells are in each cluster and then you should print it

    #add original ROI to df to check which ROIs are in each cluster
    dfFilt['file'] = roiList

    #need to add original indices to this df to check which ROIs are in each cluster; will also need to use ROI ID to pair with cell index (same index could be had by two cells from diff ROIs)
    dfFilt['index'] = idxList #idxList stores the row value of filt_df.iloc[row,column] command

    #save df to a csv
    dfFilt.to_csv(cwd+'/dfCreated/dfNeighClustered'+file[21:]+'k'+str(k)+'.csv')



def clusterCountPerROI(name):
    
    '''
    This function takes in a clustered csv and creates two csvs with the counts of each seed cell per ROI.
    One csv includes counts across all ROIs. One csv includes counts averaged across ROIs.
    Each csv has both raw count and percentage counts for each cluster.
    
    Input parameters:
        name = name of file containing clustered neighborhood data to count cluster abundance from, excluding the '.csv' (ex. 'dfNeighClusteredN60k7')
    Outputs:
        Saves two csvs saved to 'dfCreated' folder with the naming convention:
        dfClustCountsH70allk5_all.csv 
        dfClustCountsH70allk5_avg.csv
    '''
    
    import pandas as pd
    import os
    
    #get current directory
    cwd = os.getcwd()

    #read clustering csv to analyze (eg. dfNeighClusteredH70allk5; it's a csv that has each seed cell clustered)
    df = pd.read_csv(cwd+'/dfCreated/'+name+'.csv', index_col=0)
    
    dictCounts = {} #empty dict to store raw counts for each cluster
    
    #emtpy list to store patient csvs for the new df index
    dictCounts['ptList'] = [] 
    
    #create empty raw count lists; one per cluster
    for i in range(int(name[-1])): #get the number of clusters by looking at the last character in the file name
        listName = 'clust'+str(i)+'RawList'
        dictCounts[listName] = [] #give the value for each key an empty list to be appended to later
    
    #create empty percentage count lists; one per cluster
    for i in range(int(name[-1])): #get the number of clusters by looking at the last character in the file name
        listName = 'clust'+str(i)+'PercList'
        dictCounts[listName] = [] #give the value for each key an empty list to be appended to later    
    
    #for each roi in the big df
    for roi in df['file'].unique():
    
        #add roi to ptList in the correct order that matches addition of other counts to their lists
        dictCounts['ptList'].append(roi)
        
        dfSingleROI = df[df['file'] == roi] #df contains cells from only one region
    
        #raw frequency counts of cluster
        freq = dfSingleROI['cluster'].value_counts()
        
        #perc counts of cluster
        perc = dfSingleROI['cluster'].value_counts(normalize=True)
    
        #loop through each raw count in freq and each percentage count in perc
        for i in range(int(name[-1])):
            keyRaw = 'clust'+str(i)+'RawList' #create key name for each cluster; raw
            keyPerc = 'clust'+str(i)+'PercList' #create key name for each cluster; perc
            
            #if the index exists aka if the cluster has a count, add the raw and % counts
            if i in freq.index:
                dictCounts[keyRaw].append(freq[i])
                dictCounts[keyPerc].append(perc[i])
                
            #if the cluster does not have a count, add zero
            else:
                dictCounts[keyRaw].append(0)
                dictCounts[keyPerc].append(0)

    #get clinical df
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv', index_col=0)
    
    #empty list to store clinical data from dfClin
    dtrList = []
    
    #create a new df to store clinical data and scores
    dfClustCounts = pd.DataFrame()
    dfClustCounts['patient'] = dictCounts['ptList']
    dfClustCounts.set_index('patient',inplace=True) #set index labels to be the patient column
    
    #loop through each csv in ptList to get clinical values of only the csvs included in the list
    for pt in dictCounts['ptList']:
        idx = dfClin.loc[dfClin['sample'] == pt].index[0]
        dtrList.append(dfClin['dtr'].iloc[idx]) #add to dtr list    
    
    #put lists into dfScore dataframe and convert all values to int and any missing values to NaN values
    dfClustCounts['dtr'] = dtrList
    dfClustCounts['dtr'] = pd.to_numeric(dfClustCounts['dtr'],errors='coerce')   
    
    #add all scores to the df
    for key,value in dictCounts.items():
        if key != 'ptList': #don't do anything if ptList is the key because that's already been taken care of
            dfClustCounts[key] = value #add other count lists to the df
    
    #save dfClustCounts to csv
    dfClustCounts.to_csv(cwd+'/dfCreated/dfClustCounts'+name[16:]+'_all.csv')
    
    #create new dataframe with averaged numbers for every region within one patient
    dfClustCounts.index = dfClustCounts.index.str.slice(0,6)
    
    #groupby same patient (diff ROI) and take the average of the corresponding rows
    dfAvg = dfClustCounts.groupby(['patient']).mean()
    
    #save dfAvg to csv
    dfAvg.to_csv(cwd+'/dfCreated/dfClustCounts'+name[16:]+'_avg.csv')



def fig1():

    '''
    Heterogeneity across patients and tumor regions.
    This function creates figures 1b-1g. 
    All subplots are saved to the 'figures' folder.
    
    Input parameters:
        None
        
    Outputs:
        Saves csvs to 'dfCreated' folder for files containing cell densities used to generate figures.
        Saves plots to 'figures' folder for figures 1b-1g
    '''
    
    #import packages
    import plotly.graph_objs as go
    import plotly.io as pio
    import pandas as pd
    import numpy as np
    import os
    from scipy.stats import entropy,f_oneway
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.preprocessing import FunctionTransformer
    from sklearn.decomposition import PCA
    from statsmodels.stats.multicomp import pairwise_tukeyhsd


    #get dictionary to rename cells - used for subfigures
    nameDict = cellClassRename()
    
    #get current directory - used to save subfigures
    cwd = os.getcwd()
    
    #predefined color scheme; one distinct color per 11 cell types
    cellColors = {'A':'rgb(127,60,141)','B':'rgb(17,165,121)','K':'rgb(57,105,172)','C':'rgb(242,183,1)','D':'rgb(231,63,116)','J':'rgb(128,186,90)','X':'rgb(230,131,16)','E':'rgb(0,134,149)','N':'rgb(207,28,144)','H':'rgb(165,170,153)'}

    #dictionary to rename x axis labels
    xDictRename = {'pre':'P','pos':'R'}
   
    
    ##FIGURE 1b - HETEROGENEITY ACROSS PATIENTS AND TUMOR REGIONS - ALL CELL TYPES

    print('\nFigure 1b')
    print('Calculating cell densities...')

    #Calculate counts of ALL cells
    cellList = ['A','B', 'K', 'C', 'D', 'J', 'X', 'E', 'N', 'H'] #list of cells to count - exclude 'other non-immune' and 'noise' classes
    densityToCsv(cellList)
    
    print("Cell densities calculated. Csvs saved to 'dfCreated' folder.")

    #read cell density file
    df = pd.read_csv(cwd+'/dfCreated/dfDensity_all.csv',index_col=0)
    
    #empty list to store x axis labels
    xAxis = []
    
    #rename xAxis variables to 1P, 1R, etc    
    for x in df.index:
        xLabel = x[2]+xDictRename[x[3:6]]
        xAxis.append(xLabel) #add label to list
    
    traces = []
    
    for col in list(df.columns)[1:]:
        trace = go.Box(x=xAxis,y=df[col],name=nameDict[col],boxpoints='all',marker={'color':cellColors[col]},line={'color':cellColors[col]})
        traces.append(trace)
    
    layout = {'xaxis':{'title':'Tumor'},'yaxis':{'title':'Density (cells/mm^2)'},'boxmode':'group','plot_bgcolor':'rgba(0,0,0,0)'}
    fig = {'data':traces,'layout':layout}    
    pio.write_image(fig,cwd+'/figures/figure1b.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 1b saved to 'figures' folder.")    
    

    ##FIGURE 1c - HETEROGENEITY ACROSS PATIENTS AND TUMOR REGIONS - IMMUNE CELLS ONLY
    
    print('\nFigure 1c')

    #read cell density file
    df = pd.read_csv(cwd+'/dfCreated/dfDensity_all.csv',index_col=0)
        
    #only want immune cells
    df = df.drop(['N','H'],axis=1)

    #empty list to store x axis labels
    xAxis = []
    
    #rename xAxis variables to 1P, 1R, etc    
    for x in df.index:
        xLabel = x[2]+xDictRename[x[3:6]]
        xAxis.append(xLabel) #add label to list
    
    traces = []
    
    for col in list(df.columns)[1:]:
        trace = go.Box(x=xAxis,y=df[col],name=nameDict[col],boxpoints='all',marker={'color':cellColors[col]},line={'color':cellColors[col]})
        traces.append(trace)
    
    layout = {'xaxis':{'title':'Tumor'},'yaxis':{'title':'Density (cells/mm^2)'},'boxmode':'group','plot_bgcolor':'rgba(0,0,0,0)'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure1c.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 1c saved to 'figures' folder.")
    
    
    ##FIGURE 1C - KL DIVERGENCE
    
    print('\nFigure 1d')
        
    #get df of cell counts - using the perc composition with ALL cells (immune, tumor, aSMA)
    dfAll = pd.read_csv(cwd+'/dfCreated/dfDensity_all.csv',index_col=0)
    dfAvg = pd.read_csv(cwd+'/dfCreated/dfDensity_avg.csv',index_col=0)
    
    #get average of dfAvg - average values per patient (P+R)/2
    dfAvgCopy = dfAvg.copy()
    dfAvgCopy.index = dfAvgCopy.index.str.slice(2,3)
    dfAvgAvg = dfAvgCopy.groupby(['patient']).mean()
    
    #get average of primary only and recurrent only 
    dfPrim = dfAll[['pre' in i for i in dfAll.index]]
    dfPrimAvg = dfPrim.mean()
    
    dfRec = dfAll[['pos' in i for i in dfAll.index]]
    dfRecAvg = dfRec.mean()
    
    #get average of cohort
    dfCohortAvg = dfAvgAvg.mean()
    
    #get clinical df with anatomical location to add a column for site on dfAll
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv',index_col=0)
    dfClin.set_index('sample',inplace=True)
    
    dfAllCopy = dfAll.copy()
    dfAllCopy['anatomy'] = dfClin['anatomy']
    anatDict = dict(zip(dfAllCopy.index, dfAllCopy.anatomy)) #store roi:anatomy pairs for later use
    
    #get average of each anatomy region (OC,OP,L)
    dfOCAvg = dfAllCopy[dfAllCopy['anatomy'] == 'oral cavity'].mean()
    dfOPAvg = dfAllCopy[dfAllCopy['anatomy'] == 'oropharynx'].mean()
    dfLAvg = dfAllCopy[dfAllCopy['anatomy'] == 'larynx'].mean()
    
    #empty list to store KL divergence
    entListA = [] #for intra-tumoral, compares to single tumor's average
    entListB = [] #for intra-patient, compares to single patient's average (P+R)/2
    entListC = [] #for inter-patient, compares to all primary, or all recurrent depending on which one the region is from
    entListD = [] #for inter-patient, compares to cohort average across ALL pts and tumors, regardless of P or R
    entListY = [] #for intra-site, compares to same anatomic site's average across the cohort

    #create list of rois (patient samples)
    ptList = list(dfAll.index)
    
    #loop through ROIs
    for pt in ptList:
        #get row values for A-H counts for each ROI (per the pt variable)
        countRoi = dfAll.loc[pt,'A':'H'] 
        
        #get character(s) to match other df's index column
        patA = pt[0:6] #same tumor
        patB = pt[2] #same pt
        
        #get row values for A-H counts for each tumor (averaged, use the pat variables)
        countAvgA = dfAvg.loc[patA,'A':'H']
        countAvgB = dfAvgAvg.loc[patB,'A':'H']
        
        if 'pre' in pt:
            countAvgC = dfPrimAvg['A':'H']
            
        elif 'pos' in pt:
            countAvgC = dfRecAvg['A':'H']
    
        #values for matching anatomic site
        if anatDict[pt] == 'oral cavity':
            countAvgY = dfOCAvg['A':'H']
        elif anatDict[pt] == 'oropharynx':
            countAvgY = dfOPAvg['A':'H']
        elif anatDict[pt] == 'larynx':
            countAvgY = dfLAvg['A':'H']
    
        #cohort avg
        countAvgD = dfCohortAvg['A':'H']
        
        #calculate the KL Divergence using the entropy function (log2) for each avg distribution; sample distribution's divergence from the tumor's, patient's, cohort's average (qk) distribution
        entA = entropy(countRoi,qk=countAvgA,base=2) #intra-tumor
        entB = entropy(countRoi,qk=countAvgB,base=2) #intra-patient
        entC = entropy(countRoi,qk=countAvgC,base=2) #inter-patient (p or r)
        entD = entropy(countRoi,qk=countAvgD,base=2) #inter-patient (all)
        entY = entropy(countRoi,qk=countAvgY,base=2) #intra-anatomic site

        #store KL divergence in list
        entListA.append(entA)
        entListB.append(entB)
        entListC.append(entC)
        entListD.append(entD)
        entListY.append(entY)

    #create new df with KL divergence info stored in one column
    dfEnt = pd.DataFrame({'Intra-Tumor (P or R only)':entListA,'Intra-Patient (P and R)':entListB,'Inter-Patient (P or R only)':entListC,'Inter-Patient (Same Anatomic Site)':entListY,'Inter-Patient (all)':entListD},index=ptList)
    
    #empty list to store x axis labels
    xAxis = []
    #rename xAxis variables to 1P, 1R, etc
    for x in dfEnt.index:
        xLabel = x[2]+xDictRename[x[3:6]]
        xAxis.append(xLabel) #add label to list
    
    traces = []
    for col in list(dfEnt.columns):        
        trace = go.Box(x=xAxis,y=dfEnt[col],boxpoints='all',name=col)
        traces.append(trace)
    layout = {'xaxis':{'title':'Tumor'},'yaxis':{'title':'KL Divergence'},'boxmode':'group'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure1d.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 1d saved to 'figures' folder.")

    #run ANOVA to see if means are the same or different; p is significant, so not all group means are equal! 
    anP = f_oneway(entListA,entListB,entListC,entListY,entListD)[1]
    print('One way ANOVA p-value: ',anP) #one way anova

    #if ANOVA is significant, run Tukey post hoc analysis to determine which pairs are significant
    if anP < 0.05:        
        dfTukey = pd.DataFrame({'value':entListA+entListB+entListC+entListY+entListD,'label':np.repeat(list(dfEnt.columns),repeats=len(dfEnt))})
        tukey = pairwise_tukeyhsd(endog=dfTukey['value'],groups=dfTukey['label'],alpha=0.05)
        print('Tukey HSD Post-Hoc Test: ',tukey)

   
    ##FIGURE 1e - HIERARCHICAL CLUSTERING OF TUMOR REGION COMPOSITION
    
    print('\nFigure 1e')
    
    #read df
    dfFilt = pd.read_csv(cwd+'/dfCreated/dfDensity_all.csv',index_col=0)

    #get clinical csv with anatomic site
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv',index_col=0)
    dfClin.set_index('sample',inplace=True)
    
    ##color coding - three columns
    #add a new column for one color - patient number (1-9)
    dfFilt['patient'] = dfFilt.index.str.slice(2,3)
    
    #add a second new column for second color - P(pre)/R(post) (2)
    dfFilt['color2'] = dfFilt.index.str.slice(3,6)
    
    #add third column for anatomical site (3)
    dfFilt['anatomy'] = dfClin.anatomy

    #make color vector dictionary; two in total for two color palettes
    palette1 = dict(zip(dfFilt.patient.unique(),sns.color_palette('tab10',n_colors = 9))) ##9 colors for 9 patients
    palette2 = dict(zip(dfFilt.color2.unique(),sns.color_palette("Set3",len(dfFilt.index)))) #Set3 gives light green/yellow combo - adjust in illustrator
    palette3 = dict(zip(dfFilt.anatomy.unique(),sns.color_palette("Dark2",len(dfFilt.index)))) #Dark2 gives green/orange combo; Set3 gives light green/yellow combo

    #make groupings
    grouping1 = dfFilt.patient #group off of color1 column; renamed to patient so yaxis title says patient rather than color1
    grouping2 = dfFilt.color2 #group off of color2 column
    grouping3 = dfFilt.anatomy #group off anatomy column

    #map palettes onto the corresponding grouping series; names will be used as labels for color columns
    ptColors1 = pd.Series(grouping1,name='Patient ID').map(palette1)
    ptColors2 = pd.Series(grouping2,name='Primary or Recurrent').map(palette2)
    ptColors3 = pd.Series(grouping3,name='Anatomic Site').map(palette3)
    dfColors = pd.concat([ptColors1,ptColors2,ptColors3], axis=1)    
    
    #get only the desired columns aka drop dtr, added three cols for palettes
    dfFiltSubset = dfFilt[list(dfFilt.columns)[1:-3]]
    
    #rename columns from cell classes to words
    columns = dfFiltSubset.columns #get columns to rename
    #update columns to use real cell classes rather than the IDs
    newCols = []
    for c in columns:
        newCols.append(nameDict[c])
    dfFiltSubset.columns = newCols
    
    #log10+1 scale
    dfFiltLog = np.log10((dfFiltSubset)+1)
    
    #plot the clusttermap with the colors
    graph = sns.clustermap(dfFiltLog,method='ward',metric='euclidean',cmap='vlag',row_colors=dfColors,yticklabels=True) #yticklabels = true to show all patient labels; row_colors are values of ptColors
    
    #add x-axis label and title to graph
    ax = graph.ax_heatmap
    ax.set_xlabel("\nCell Phenotype")
    ax.set_ylabel("Tumor Region")
    
    #adjust y limit to fix bug in matplotlib code: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    bottom,top = ax.get_ylim()
    ax.set_ylim(bottom+.5,top-.5)
    
    plt.savefig(cwd+'/figures/figure1e.png',format='png',bbox_inches='tight') #saves figure to 'figures' folder as a png file
    plt.close()
    
    print("Figure 1e saved to 'figures' folder.")
    
    
    ##FIGURE 1f - PRINCIPAL COMPONENT ANALYSIS BY PATIENT/TIMEPOINT
    
    print('\nFigure 1f')
    
    #read dataframe
    df = pd.read_csv(cwd+'/dfCreated/dfDensity_all.csv',index_col=0)
    
    #set color palette so it matches other plots; patient colors
    palette = dict(zip(df.index.str.slice(2,3).unique(),sns.color_palette('tab10',n_colors = 9))) ##9 colors for 9 patients
    
    #get roiList, for later dataframe creation, to map to
    roiList = list(df.index)
        
    #separate out features; drop dtr column
    feat = df.iloc[:,1:].values
    
    #standardize the features w log10+1
    transformer = FunctionTransformer(np.log10)
    feat = transformer.fit_transform((feat)+1) 
    
    #target = patient ID, P/R
    df['target'] = [roi[2] for roi in df.index.values]
    df['time'] = [roi[3:6] for roi in df.index.values] #to plot P vs R as different style markers (style param in scatter line)
    
    #do PCA
    pca = PCA(n_components = 2) #generate model based on top two PCs
    principalComponents = pca.fit_transform(feat) #fit and transform features
    
    #create new dataframe holding the principal components; need to specify index so that you can map the target column in next step
    dfPrinc = pd.DataFrame(data=principalComponents, columns=['PC1','PC2'], index = roiList)
    
    #concatenate the principal component df and the target column (this is what gets color coded in final plot)
    dfFinal = pd.concat([dfPrinc,df[['target','time']]], axis=1)
    
    #rename target column name to appropriate name for plotting - to patient
    dfFinal = dfFinal.rename(columns={'target':'Patient'})
    
    #visualize using seaborn, color according to the target column
    sns.set(rc={'figure.figsize':(6,6)})
    graph = sns.scatterplot(dfFinal['PC1'], dfFinal['PC2'], hue=dfFinal['Patient'], legend='full',palette=palette,style=dfFinal['time'])
    graph.legend(loc='lower right',bbox_to_anchor=(1.25, 0.5),ncol=1)
    
    plt.savefig(cwd+'/figures/figure1f.png',format='png',bbox_inches='tight') #saves figure to 'figures' folder as a png file
    plt.close()    
    
    print("Figure 1f saved to 'figures' folder.")
        
    
    ##FIGURE 1g - PRINCIPAL COMPONENT ANALYSIS BY ANATOMIC SITE
    
    print('\nFigure 1g')
    
    #read dataframe
    df = pd.read_csv(cwd+'/dfCreated/dfDensity_all.csv',index_col=0)
    
    #get clinical df with anatomical location
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv',index_col=0)
    dfClin.set_index('sample',inplace=True)    
        
    #get roiList, for later dataframe creation, to map to
    roiList = list(df.index)
        
    #separate out features; drop dtr column
    feat = df.iloc[:,1:].values
    
    #standardize the features w log10+1
    transformer = FunctionTransformer(np.log10)
    feat = transformer.fit_transform((feat)+1) 

    #target = anatomic site
    df['target'] = dfClin['anatomy']
    
    #set color palette so it matches other plots; patient colors
    palette = dict(zip(df.target.unique(),sns.color_palette('tab10',n_colors = 3))) #3 colors for 3 sites
    
    #do PCA
    pca = PCA(n_components = 2) #generate model based on top two PCs
    principalComponents = pca.fit_transform(feat) #fit and transform features
    
    #create new dataframe holding the principal components; need to specify index so that you can map the target column in next step
    dfPrinc = pd.DataFrame(data=principalComponents, columns=['PC1','PC2'], index = roiList)
    
    #concatenate the principal component df and the target column (this is what gets color coded in final plot)
    dfFinal = pd.concat([dfPrinc,df[['target']]], axis=1)
    
    #rename target column name to appropriate name for plotting - to site
    dfFinal = dfFinal.rename(columns={'target':'Anatomic Site'})
    
    #visualize using seaborn, color according to the target column
    sns.set(rc={'figure.figsize':(6,6)})
    graph = sns.scatterplot(dfFinal['PC1'], dfFinal['PC2'], hue=dfFinal['Anatomic Site'], legend='full',palette=palette)
    graph.legend(loc='lower right',bbox_to_anchor=(1.25, 0.5),ncol=1)
    
    plt.savefig(cwd+'/figures/figure1g.png',format='png',bbox_inches='tight') #saves figure to 'figures' folder as a png file
    plt.close()    
    
    print("Figure 1g saved to 'figures' folder.")
  


def fig2():
    
    '''
    Tumor cellular composition changes following therapy.
    This function creates figures 2a-2d. 
    All subplots are saved to the 'figures' folder.
    
    Note: This function assumes densityToCsv() has already been run (this happens in fig1() function).
    
    Input parameters:
        None
        
    Outputs:
        Saves plots to 'figures' folder for figures 2a-2d
    '''
    
    #import packages
    import pandas as pd
    import os
    import plotly.graph_objs as go
    import plotly.io as pio
    from scipy.stats import ttest_rel, wilcoxon
    from statsmodels.stats.multitest import fdrcorrection
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.preprocessing import MinMaxScaler

    #get dictionary to rename cells - used for subfigures
    nameDict = cellClassRename()
    
    #get current directory - used to save subfigures
    cwd = os.getcwd()


    ##FIGURE 2a - PRIMARY VS RECURRENT CELL COUNTS
    
    print('\nFigure 2a')
    
    #read in counts df
    df = pd.read_csv(cwd+'/dfCreated/dfDensity_avg.csv',index_col=0)
    
    #drop the dtr column
    dfFilt = df.iloc[:,1:]
    
    #separate out pre and post treatment tumors
    dfPre = dfFilt[['pre' in row for row in dfFilt.index]]
    dfPost = dfFilt[['pos' in row for row in dfFilt.index]]
    
    #put pre and post df's in dfList
    dfList = [dfPre,dfPost]
    
    #get corresponding cell names from cellClassRename function
    nameDict = cellClassRename()
    
    #empty list to store two traces for plotting
    traces = []
    
    #loop through the two dfs
    for i in range(len(dfList)):
    
        #get df
        df = dfList[i]
    
        #empty lists to store data for plotting - one per trace
        xAxis = [] #will store label names
        yList = [] #will store count values to plot in boxplot
    
        for col in df.columns: #for each column
    
            for row in df.index: #for each row of patients
                xAxis.append(nameDict[col]) #add the column's real class name to the xAxis label list (each label should repeat itself nine times)
                yList.append(df.loc[row,col]*100) #add the column/row's value (aka actual value to plot)
    
        #naming of traces
        if i == 0:
            treat = 'Primary'
        elif i == 1:
            treat = 'Recurrent'
    
        #create one trace for pre and one for post
        trace = go.Box(x=xAxis,y=yList,name=treat,boxmean=True,boxpoints='all')
        traces.append(trace)
    
    layout = {'xaxis':{'title':'Cell Phenotype'},'yaxis':{'title':'Density (cells/mm^2) (Averaged Across Regions)'},'boxmode':'group'}    
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure2a.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Figure 2a saved to 'figures' folder.")
    
    #check for significance between pre/post tumors and print the adjusted p-value, paired tests
    print('Calculating p-values...')
    
    pList = [] #empty list to store p values for MHT correction
    
    for col in dfPre.columns: #for each column in the df
        valueList = [] #overall list to store two lists to compare (pre and post); reset with each new column
        for df in dfList: #for each df (pre and then post)
            valList = list(df[col])
            valueList.append(valList) #add pre or post's value list to overall valueList
                
        #non-parametric paired t-test (for tumor densities, whose difference does NOT follow a normal distribution)
        if col == 'H':
            pVal = wilcoxon(valueList[0],valueList[1])[1]
        
        #all other cell types' differences DO follow a normal distribution, so use parametric paired t-test
        else:
            pVal = ttest_rel(valueList[0],valueList[1])[1]
            
        pList.append(pVal)
        
    #correct for MHT
    mhtList = fdrcorrection(pList,alpha=0.05)[1]
    
    print('\nMHT corrected p-values:')
    for i in range(len(dfPre.columns)):
        col = dfPre.columns[i]
        mhtP = mhtList[i]
        print(nameDict[col],': ',mhtP)        
        

    ##FIGURE 2b - PRIMARY VS RECURRENT CELL COUNTS

    print('\nFigure 2b')

    #read in counts df
    df = pd.read_csv(cwd+'/dfCreated/dfDensity_avg.csv',index_col=0)
    
    #drop the dtr and dtd columns
    dfFilt = df.iloc[:,1:]
    
    #separate out pre and post treatment tumors
    dfPre = dfFilt[['pre' in row for row in dfFilt.index]]
    dfPost = dfFilt[['pos' in row for row in dfFilt.index]]
    
    #predefined list of patients 1-9
    colList = ['1','2','3','4','5','6','7','8','9']
    
    #make title list for each subplot
    cols = list(dfFilt.columns)
    titleList = []
    for c in cols:
        
        if c != 'N':
            x = nameDict[c]
        else:
            x = 'aSMA+ Mesenchymal' #matplot cannot print the alpha symbol so renaming for this plot to an 'a'
        titleList.append(x)
    
    #update dfPre
    dfPreT = dfPre.T
    dfPreT.columns = colList #update col names
    
    #update row names
    rowNames = list(dfPreT.index)
    rowUpdate = []
    
    for r in rowNames:
        row = nameDict[r]+'-P'
        rowUpdate.append(row)    
    dfPreT.index = rowUpdate    
    
    #update dfPost
    dfPostT = dfPost.T
    
    dfPostT.columns = colList #update col names
    
    #update row names
    rowNames = list(dfPostT.index)
    rowUpdate = []
    
    for r in rowNames:
        row = nameDict[r]+'-R'
        rowUpdate.append(row)  
    dfPostT.index = rowUpdate
     
    #merge dfs by alternating rows
    s1 = dfPreT.assign(a = np.arange(len(dfPreT))).set_index('a', append=True) 
    s2 = dfPostT.assign(a = np.arange(len(dfPostT))).set_index('a', append=True)
    
    dfFinal = (pd.concat([s1, s2], keys=('Primary','Recurrent'))
            .sort_index(kind='merge', level=2)
            .reset_index(level=2, drop=True)
            .swaplevel(0,1))
    
    #reset the index to just primary and recurrent
    dfIdx = dfFinal.index.droplevel(level=0)
    dfFinal.index = dfIdx
    
    #set up plotting
    fig, ax =plt.subplots(2,5,figsize=(25,25))
    
    #range through 11 to set variables for placement of graph, where in the df to get data from
    for n in range(11):
        
        #set placement of graph
        #row
        if n < 5:
            a = 0
        else:
            a = 1
        #column
        if n < 5:
            b = n
        else:
            b = n-5
        
        #set iloc to draw data from
        i = n*2
        j = i+1
        
        if i < 20:
            
            ax[a,b].set_title(titleList[n],fontdict={'fontsize':14}) #set title of each subplot
            ax[a,b].set_ylabel("Density (cells/mm^2) (Averaged Across Regions)") #set y axis label of each subplot        
            #make lineplot; use tab10 color scheme to match before
            sns.lineplot(data=dfFinal.iloc[[i,j],:],dashes=False,palette='tab10',ax=ax[a,b],legend=False,marker="o")
    
        #skip the final plot that does not exist 
        else:
            break
    
    #final adjustments
    ax[1,4].legend(['1', '2', '3','4','5','6','7','8','9'],bbox_to_anchor=(1,0.75),ncol=1,title='Patient',prop={'size': 18}) #create legend for the final ax 
    
    plt.savefig(cwd+'/figures/figure2b.png',format='png') #saves figure to 'figures' folder as a png file
    plt.close()
    print("Figure 2b saved to 'figures' folder.")


    ##FIGURE 2c - HIERARCHICAL CLUSTERING OF AVERAGE CHANGE IN TIME CELLULAR COMPOSITION
    
    print('\nFigure 2c')
    
    #read df
    df = pd.read_csv(cwd+'/dfCreated/dfDensity_avg.csv',index_col=0)
    
    #drop the dtr column
    dfFilt = df.iloc[:,1:]
    
    #separate out pre and post treatment tumors
    dfPre = dfFilt[['pre' in row for row in dfFilt.index]]
    dfPost = dfFilt[['pos' in row for row in dfFilt.index]]
    
    #put pre and post df's in dfList
    dfList = [dfPre,dfPost]
    
    #empty dict to store col:difference list for each cell type
    diffDict = {}
    
    #get a diffList for each column in the df; these will be the features to cluster on
    for col in dfPre.columns: #for each column in the df
        valueList = [] #overall list to store two lists to compare (pre and post); reset with each new column
        for df in dfList: #for each df (pre and then post)
            valList = list(df[col])
            valueList.append(valList) #add pre or post's value list to overall valueList
    
        diffList = [x1 - x0 for (x1, x0) in zip(valueList[1], valueList[0])]
    
        if col != 'N':
            diffDict[nameDict[col]] = diffList
        else:
            diffDict['aSMA+ Mesenchymal'] = diffList #matplot cannot print the alpha symbol so renaming for this plot to an 'a'
    
    dfDiff = pd.DataFrame(diffDict,index=['1','2','3','4','5','6','7','8','9'])
    
    #get treatment data and add new column to dfDiff for color coding
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv', index_col=0)
    dfClin = dfClin.dropna() #drop recurrent tumors
    dfClin.set_index('sample',inplace=True)
    dfClin.index = dfClin.index.str.slice(0,6) #update patient column with just the first six characters in string; so you can use groupby
    dfClin = dfClin[~dfClin.index.duplicated(keep='first')]
    dfClin.index = dfClin.index.str.slice(2,3)
    
    dfDiff['patient'] = dfDiff.index
    dfDiff['tx'] = dfClin['tx']
        
    #make color vector dictionary; just one needed for 1 color palette (patient ID)
    palette1 = dict(zip(dfDiff.patient,sns.color_palette("tab10",n_colors = 9))) #9 colors for 9 patients
    palette2 = dict(zip(dfDiff.tx.unique(),sns.color_palette("Accent",n_colors=3))) #3 colors becauuse 3 tx options
    
    #make groupings
    grouping1 = dfDiff.patient #group off of patient column
    grouping2 = dfDiff.tx #group off of tx column
    
    #map palettes onto the corresponding grouping series
    ptColors1 = pd.Series(grouping1,name='Patient ID').map(palette1)
    ptColors2 = pd.Series(grouping2,name='Therapy').map(palette2)
    dfColors = pd.concat([ptColors1,ptColors2], axis=1)
    
    #drop the last two columns before clustering on values
    dfDiff = dfDiff[list(dfDiff.columns)[:-2]]
    
    ##normalize [-1,1] while accounting for negative values
    #create an all positive values of dfDiff
    dfDiffPos = dfDiff.abs()
    
    #scale differences [0,1]
    xArray = dfDiffPos.values #return numpy array of all values in the df  
    columns = dfDiffPos.columns #get cols to add back to dfNorm
    patients = dfDiffPos.index
    
    scaler = MinMaxScaler(feature_range=(0,1))
    x_scaled = scaler.fit_transform(xArray)
    dfNorm = pd.DataFrame(x_scaled,index=patients,columns=columns)
    dfNew = pd.DataFrame(index=patients,columns = columns) #create empty df with same indices and columns to be filled w scaled data

    #fill in dfNew with appropriate +/- values from dfNorm
    for c in range(len(dfDiff.columns)):
        for i in range(len(dfDiff)):
            if dfDiff.iloc[i,c] >= 0: #if positive or zero add value as is
                dfNew.iloc[i,c] = dfNorm.iloc[i,c] 
            else: #if negative
                dfNew.iloc[i,c] = dfNorm.iloc[i,c]*(-1) #then flip the sign of the normalized value
                #some negative zeros but they equal positive zeros: https://stackoverflow.com/questions/4083401/negative-zero-in-python?lq=1
    
    #convert dfNew to floats
    dfNewFloat = dfNew.apply(pd.to_numeric)
    
    #plot the clusttermap with the colors
    graph = sns.clustermap(dfNewFloat,method='ward',metric='euclidean',cmap='vlag',row_colors=dfColors,yticklabels=True)
    
    #add x-axis label and title to graph
    ax = graph.ax_heatmap
    ax.set_xlabel("\nCell Phenotype")
    ax.set_ylabel("Patient")
    
    #adjust y limit to fix bug in matplotlib code
    bottom,top = ax.get_ylim()
    ax.set_ylim(bottom+.5,top-.5)
    
    plt.savefig(cwd+'/figures/figure2c.png',format='png',bbox_inches='tight') #saves figure to 'figures' folder as a png file
    #plt.show(graph)
    plt.close()
    
    print("Figure 2c saved to 'figures' folder.")
                

    ##FIGURE 2d - CHANGE IN CELL DENSITIES WITH PFS COLOR CODED
    
    print('\nFigure 2d')

    survColor = {'Short':'rgb(235,173,5)','Long':'rgb(123,60,171)'}
    
    #load df
    df = pd.read_csv(cwd+'/dfCreated/dfDensity_avg.csv',index_col=0)
    
    #drop the dtr column
    dfFilt = df.iloc[:,1:]
    
    #separate out pre and post treatment tumors
    dfPre = dfFilt[['pre' in row for row in dfFilt.index]]
    dfPost = dfFilt[['pos' in row for row in dfFilt.index]]
    
    #put pre and post df's in dfList
    dfList = [dfPre,dfPost]
    
    #empty dict to store col:difference list for each cell type
    diffDict = {}
    
    #get a diffList for each column in the df; these will be the features to cluster on
    for col in dfPre.columns: #for each column in the df
        valueList = [] #overall list to store two lists to compare (pre and post); reset with each new column
        for df in dfList: #for each df (pre and then post)
            valList = list(df[col])
            valueList.append(valList) #add pre or post's value list to overall valueList
    
        diffList = [x1 - x0 for (x1, x0) in zip(valueList[1], valueList[0])]
    
        diffDict[nameDict[col[0]]] = diffList
    
    dfDiff = pd.DataFrame(diffDict,index=['1','2','3','4','5','6','7','8','9'])
    dfDiff.index.name = 'patient'
    
    #get clinical df and add dtr to dfDiff
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv', index_col=0)
    dfClin.set_index('sample',inplace=True) #set index labels to be the patient column
    dfClinPre = dfClin[['pre' in row for row in dfClin.index]]
    
    dfClinPre.index = dfClinPre.index.str.slice(2,3)
    dfClinPreCopy = dfClinPre.copy()
    dfClinPreCopy['dtr'] = pd.to_numeric(dfClinPreCopy['dtr'],errors='coerce')
    
    dfClinAvg = dfClinPreCopy.groupby(['sample']).mean()
    dfDiff['dtr'] = dfClinAvg['dtr']
    
    #split pts by median survival
    med = dfDiff.median()['dtr']
    
    dfShort = dfDiff[dfDiff['dtr'] < med]
    dfLong = dfDiff[dfDiff['dtr'] >= med]
    dfList = [dfShort,dfLong]
    
    traces = []
    for i in range(len(dfList)):
        #get short or long term df
        df = dfList[i]
        
        #empty lists to store data for plotting - one per trace
        xAxis = [] #will store label names
        yList = [] #will store delta values to plot in boxplot
    
        for col in df.columns[:-1]: #for each column, drop dtr
            for row in df.index: #for each row of patients
                xAxis.append(col) #add the column's real class name to the xAxis label list (each label should repeat itself nine times)
                yList.append(df.loc[row,col]) #add the column/row's value (aka actual value to plot)
                    
        #name trace
        if i == 0:
            survival = 'Short'
        elif i ==1:
            survival = 'Long'
    
        #create one trace for short/long
        trace = go.Box(x=xAxis,y=yList,name=survival+' PFS',boxmean=True,boxpoints='all',marker={'color':survColor[survival]},fillcolor='rgba(255,255,255,0)',line={'color':'rgba(255,255,255,0)'})
        traces.append(trace)
        
    #reset xAxis and yList 
    xAxis = []
    yList = []
    #a separate box plot of deltas for each cell type for the box, no color for the points
    for col in dfDiff.columns[:-1]: #for each column, drop dtr
        for row in dfDiff.index: #for each row of patients
            xAxis.append(col)
            yList.append(dfDiff.loc[row,col])
    
    traceBox = go.Box(x=xAxis,y=yList,boxmean=True,boxpoints=False,fillcolor='rgb(232,232,232)',line={'color':'rgb(178,178,178)'},name='Change')
    traces.append(traceBox)
    layout = {'xaxis':{'title':'Cell Phenotype'},'yaxis':{'title':'Change in Density (cells/mm^2) (Averaged Across Regions)'}}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure2d.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Figure 2d saved to 'figures' folder.")



def fig3():
    
    '''
    Mixing score quantifies the spatial organization of tumors.
    This function creates figures 3b-3i. 
    All subplots are saved to the 'figures' folder.
    
    Note: This function assumes densityToCsv() has already been run (this happens in fig1() function).
    
    Input parameters:
        None
        
    Outputs:
        Saves csvs to 'dfCreated' folder for files containing mixing scores, PD-1+ counts and Ki-67+ counts used to generate figures.
        Saves plots to 'figures' folder for figures 3b-3i
    '''
  
    #import packages    
    import pandas as pd
    import numpy as np
    import os
    import plotly.graph_objs as go
    import plotly.io as pio
    from scipy.stats import mannwhitneyu, ttest_ind,f_oneway
    from sklearn.preprocessing import FunctionTransformer
    from sklearn.decomposition import PCA
    import seaborn as sns
    import matplotlib.pyplot as plt
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
    from statsmodels.stats.multitest import fdrcorrection

    #get dictionary to rename cells - used for subfigures
    nameDict = cellClassRename()
    
    #get current directory - used to save subfigures
    cwd = os.getcwd()
    
    #set color dict to map colors from for plotting spatial organization
    colorDict = {'Cold':'#E377C2','Mixed':'#FF7F03','Compartmentalized':'#2CA02C'}

    #rename xAxis variables to 1P, 1R, etc
    xDictRename = {'pre':'P','pos':'R'}

    #create function to map spatial labels into dfs within figure 3 subplots
    def map_mix(mix):
    
        if mix < 0.107: #value corresponds to median of primary tumor regions
            return 'Compartmentalized'
        elif mix > 0.107:
            return 'Mixed'

    
    ##FIGURE 3b - MIXING SCORE ACROSS TUMOR REGIONS
    
    print('\nFigure 3b')
    
    #calculate mixing score for all tumor regions and save data to new csv
    print('Calculating mixing scores...')
    calculateMixing(30,['A','B','C','D','E','J','K','X']) #30px = 15µm radius, all immune cells
    print("Mixing scores calculated. Csvs saved to 'dfCreated' folder.")

    #determine which tumor regions are cold
    print('Calculating which tumor regions are cold...')
    calculateCold()
    print("New mixing score csvs excluding cold regions saved to 'dfCreated' folder.")
    
    #read csv
    df = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    
    #use map_mix function defined above to map the appropriate label to the 'spatiial' column of df
    df['spatial'] = df['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold ROIs
    coldList = ['pt2post_roi1','pt9post_roi3'] #this is a pre-defined list based on # of immune cells from calculateCold()
    for roi in coldList:
            df.loc[roi,'spatial'] = 'Cold'
    
    #rename xaxis labels to be primary recurrent rather than pre/post
    ptList = []
    for x in df.index:
        xLabel = x[2]+xDictRename[x[3:6]]
        ptList.append(xLabel) #add label to list
    df.index = ptList
    
    #get tumor list
    xAxis = list(df.index)
    
    #trace with all data points, box only
    traceData = go.Box(x=xAxis,y=df['mixScore'],boxmean=True,boxpoints=False,name='data',jitter=1)
    
    #create traces to store all trace variables
    traces = [traceData]
    
    #create new trace of markers only for color-coding by spatial org
    for cat in df['spatial'].unique():
    
        #get a sub-df with rois that all have the same spatial category
        dfSpat = df[df['spatial'] == cat]
    
        #create marker trace for that spatial category and add it to traces list
        traceSpat = go.Box(x=list(dfSpat.index),y=dfSpat['mixScore'],boxpoints='all',fillcolor='rgba(255,255,255,0)',line={'color':'rgba(255,255,255,0)'},marker={'color':colorDict[cat]},jitter=1)
        traces.append(traceSpat)
    
    #create annotations with overall averaged mixing score label
    #get overall mixing score designation from average csv 
    file2 = 'dfMixingScoresExcludeCold_avg'
    dfAvg = pd.read_csv(cwd+'/dfCreated/'+file2+'.csv', index_col=0)
    #add spatial column
    dfAvg['spatial'] = dfAvg['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold tumors where the regions are all cold
    coldList = ['pt2post'] #this is a pre-defined list based on # of immune cells from calculateCold()
    #add row with cold tumor
    for roi in coldList:
            dfAvg.loc[roi,'spatial'] = 'Cold'
    
    #update index values to be 1P, 1R, etc
    ptList = []
    for x in dfAvg.index:
        xLabel = x[2]+xDictRename[x[3:6]]
        ptList.append(xLabel) #add label to list
    dfAvg.index = ptList
    dfAvg = dfAvg.sort_index()
            
    #get lists of df columns 
    ptList = list(dfAvg.index)
    spatList = list(dfAvg['spatial'])
    mixList = list(dfAvg['mixScore'])
    
    #empty list to store tuples
    annotTupList = []
    
    #create tuples to add to the annotTupList
    for i in range(len(ptList)):
        
        if mixList[i] > 0:
            annot = (i,mixList[i]+1.5,spatList[i])
        
        else: #for nan values
            annot = (i,1.5,spatList[i])
        annotTupList.append(annot) 
    
    #empty list to store dictionaries for annotations on plot
    annotations = []
    
    #create annotation list
    for j in range(len(annotTupList)):
        annotDict = dict(x=annotTupList[j][0],y=annotTupList[j][1],text=annotTupList[j][2],showarrow=False,textangle=-90,font={'color':colorDict[annotTupList[j][2]]})
        annotations.append(annotDict)
        
    #figure layout
    layout = {'xaxis':{'title':'Tumor'},'yaxis':{'title':'Mixing Score'},'annotations':annotations}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure3b.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Figure 3b saved to 'figures' folder.")
    
    
    ##FIGURE 3c - MIXING SCORE BY ANATOMIC SITE
    
    print('\nFigure 3c')
    
    #read csv with mixing scores
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv',index_col=0)
    
    #add spatial designation
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold ROIs - all ROIs
    coldList = ['pt2post_roi1','pt9post_roi3'] #this is a pre-defined list based on # of immune cells/area
    
    for roi in coldList:
            dfMix.loc[roi,'spatial'] = 'Cold'
    
    #get clinical df with anatomical location
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv',index_col=0)
    dfClin.set_index('sample',inplace=True)
    dfMix['anatomy'] = dfClin['anatomy']
    
    #separate out regions based on anatomic site
    dfOC = dfMix[['oral cavity' in row for row in dfMix['anatomy']]]
    dfO = dfMix[['oropharynx' in row for row in dfMix['anatomy']]]
    dfL = dfMix[['larynx' in row for row in dfMix['anatomy']]]
    
    traceOC = go.Box(y=dfOC.loc[:,'mixScore'],name='Oral Cavity',boxpoints='all')
    traceO = go.Box(y=dfO.loc[:,'mixScore'],name='Oropharynx',boxpoints='all')
    traceL = go.Box(y=dfL.loc[:,'mixScore'],name='Larynx',boxpoints='all')
    traces = [traceOC,traceO,traceL]
    layout = {'xaxis':{'title':'Anatomic Site'},'yaxis':{'title':'Mixing Score'},'boxmode':'group'}
    fig = {'data':traces,'layout':layout}
    
    pio.write_image(fig,cwd+'/figures/figure3c.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 3c saved to 'figures' folder.")

    #run ANOVA to see if means are the same or different
    anP = f_oneway(dfOC.loc[:,'mixScore'],dfO.loc[:,'mixScore'],dfL.loc[:,'mixScore'])[1]
    print('One way ANOVA p-value: ',anP) #one way anova

    #if ANOVA is significant, run Tukey post hoc analysis to determine which pairs are significant
    if anP < 0.05:        
        dfTukey = pd.DataFrame({'value':dfMix['mixScore'],'label':dfMix['anatomy']})        
        tukey = pairwise_tukeyhsd(endog=dfTukey['value'],groups=dfTukey['label'],alpha=0.05)
        print('Tukey HSD Post-Hoc Test: ',tukey)


    ##FIGURE 3d - SURVIVAL CURVE, MIXING SCORE
    
    print('\nFigure 3d') 
    
    #get df from saved csv; dependent on the path
    df = pd.read_csv(cwd+'/dfCreated/dfMixingScoresExcludeCold_avg.csv', index_col=0)
    
    #filter df so only primary tumors are analyzed
    dfFilt = df[['pre' in row for row in df.index]]
      
    #remove rows with NaN values in the measure being plotted according to the survType also being plotted
    dfFilt = dfFilt.dropna(subset=['mixScore','dtr'])
    
    #median of the desired score
    med = dfFilt.median()['mixScore']
    
    #split into low and high scores
    dfLow = dfFilt[dfFilt['mixScore'] < med] #low group
    dfHigh = dfFilt[dfFilt['mixScore'] >= med] #high group
    
    #fit KM model
    kmf = KaplanMeierFitter()
    
    #create subplot to hold both curves
    plt.figure() #new figure for each call to this function
    ax = plt.subplot(111)
    
    #dfLow curve
    kmf.fit(dfLow['dtr'], event_observed=[1]*len(dfLow), label="Compartmentalized")
    kmf.plot(ax=ax,ci_show=False,c='#2CA02C') #set color to match other figures
    
    #dfHigh curve
    kmf.fit(dfHigh['dtr'], event_observed=[1]*len(dfHigh), label="Mixed")
    kmf.plot(ax=ax,ci_show=False,c='#FF7F03') #set color to match other figures
    
    plt.ylim(0, 1) #set y axis to go from zero to one
    plt.margins(0) #make (0,0) start in bottom left corner rather than a little bit off
    
    #run logrank test for p-value calculation
    results = logrank_test(dfLow['dtr'],dfHigh['dtr'],[1]*len(dfLow),[1]*len(dfHigh),alpha=0.95)
    
    plt.ylabel('Progression Free Survival (%)')
    plt.xlabel('Days')
    plt.text(x=1150,y=0.7,s='p = '+str(round(results.p_value,3)))

    plt.savefig(cwd+'/figures/figure3d.png',format='png') #saves figure to 'figures' folder as a png file
    plt.close()
    
    print("Figure 3d saved to 'figures' folder.")        


    ##FIGURE 3e - PRINCIPAL COMPONENT ANALYSIS BY MIXING SCORE
    
    print('\nFigure 3e')
    
    #set color palette so it matches other plots; mixing score colors
    palette = dict(Compartmentalized = '#2CA02C',Mixed = '#FF7F03',Cold = '#E377C2')
    
    #read counts csv
    df = pd.read_csv(cwd+'/dfCreated/dfDensity_all.csv',index_col=0)
    
    #get roiList, for later dataframe creation, to map to
    roiList = list(df.index)
    
    #separate out features; drop dtr column
    feat = df.iloc[:,1:].values
    
    #standardize the features
    transformer = FunctionTransformer(np.log10)
    feat = transformer.fit_transform((feat)+1)
        
    #add a column on mixing df for mixed, compartmentalized, cold - to see intra-tumor heterogeneity    
    #read in mixing score csv - all regions
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix)) #map_mix defined above
    
    #adjust df to also account for cold ROIs - all ROIs
    coldList = ['pt2post_roi1','pt9post_roi3'] #this is a pre-defined list based on # of immune cells from calculateCold()
    for roi in coldList:
            dfMix.loc[roi,'spatial'] = 'Cold'
    
    #now add the spatial header into the target column in the df
    df['target'] = list(dfMix['spatial'])
    
    #do PCA
    pca = PCA(n_components = 2) #create model based on top two PCs
    principalComponents = pca.fit_transform(feat) #fit and transform features
    
    #create new dataframe holding the principal components; need to specify index so that you can map the target column in next step
    dfPrinc = pd.DataFrame(data=principalComponents, columns=['PC1','PC2'], index = roiList)
    
    #concatenate the principal component df and the target column (this is what gets color coded in final plot)
    dfFinal = pd.concat([dfPrinc,df[['target']]], axis=1)
    
    #rename target column name to appropriate name for plotting
    dfFinal = dfFinal.rename(columns={'target':'Spatial Organization'})
    
    #visualize using seaborn, color according to the target column
    sns.set(rc={'figure.figsize':(7,7)})
    sns.scatterplot(dfFinal['PC1'], dfFinal['PC2'], hue=dfFinal['Spatial Organization'], legend='full',palette=palette)
    
    plt.savefig(cwd+'/figures/figure3e.png',format='png') #saves figure to 'figures' folder as a png file
    plt.close()
    
    print("Figure 3e saved to 'figures' folder.")    

    
    ##FIGURE 3f - CELL DENSITIES BY MIXING SCORE
    
    print('\nFigure 3f')

    #read in density df
    df = pd.read_csv(cwd+'/dfCreated/dfDensity_all.csv',index_col=0)
    
    #drop the dtr columns
    dfFilt = df.iloc[:,1:]
    
    #read in mix df
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold ROIs - all ROIs
    coldList = ['pt2post_roi1','pt9post_roi3'] #this is a pre-defined list based on # of immune cells/area
    for roi in coldList:
            dfMix.loc[roi,'spatial'] = 'Cold'
    
    #now add the spatial header into the spatial column to the counts df
    dfFilt['spatial'] = list(dfMix['spatial'])
    
    #separate out mixed and compartmentalized based on spatial column
    dfMixed = dfFilt[['Mixed' in row for row in dfFilt['spatial']]]
    dfCompart = dfFilt[['Compartmentalized' in row for row in dfFilt['spatial']]]
    
    #put M and C df's in dfList
    dfList = [dfCompart,dfMixed]
        
    #empty list to store two traces for plotting
    traces = []
    
    #loop through the two dfs
    for i in range(len(dfList)):
        #get df
        df = dfList[i]
    
        #empty lists to store data for plotting - one per trace
        xAxis = [] #will store label names
        yList = [] #will store count values to plot in boxplot
    
        for col in df.columns[:-1]: #for each column
    
            for row in df.index: #for each row of patients
                xAxis.append(nameDict[col[0]]) #add the column's real class name to the xAxis label list (each label should repeat itself nine times)
                yList.append(df.loc[row,col]) #add the column/row's value (aka actual value to plot)
    
        #naming of traces
        if i == 0:
            spat = 'Compartmentalized'
        elif i == 1:
            spat = 'Mixed'
    
        #create one trace for C and M for post
        trace = go.Box(x=xAxis,y=yList,name=spat,boxmean=True,boxpoints='all',marker={'color':colorDict[spat]},line={'color':colorDict[spat]})
        traces.append(trace)
    
    layout = {'xaxis':{'title':'Cell Phenotype'},'yaxis':{'title':'Density (cells/mm^2)'},'boxmode':'group'}
    fig = {'data':traces,'layout':layout}
    
    pio.write_image(fig,cwd+'/figures/figure3f.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    print("Figure 3f saved to 'figures' folder.")
   
    #p-value significance testing
    print('Calculating p-values...')
    
    pDict = {}
    for col in dfMixed.columns[:-1]:        
        if col in ['C','J','X,','H']: #for cells whose differences DO NOT follow a normal distribution; use non-parametric independent t-test aka M-W U test
            #adjust for different alternatives 
            if col == 'H':
                pOne = mannwhitneyu(dfCompart.loc[:,col],dfMixed.loc[:,col],alternative='less')[1] #alt = compart is LESS than mixed for tumor only
            else:
                pOne = mannwhitneyu(dfCompart.loc[:,col],dfMixed.loc[:,col],alternative='greater')[1] #alt = compart is GREATER than mixed for other cell types
            
        else: #for cells whose differences DO follow a normal distribution; use parametric independent t-test
            pTwo = ttest_ind(dfCompart.loc[:,col],dfMixed.loc[:,col],equal_var=False)[1] #returns results from TWO-TAIL independent (parametric) t-test        
            pOne = pTwo/2 #divide by two to get one tailed p value
        
        #add cell type:p values to dict
        pDict[col] = pOne
        
    ##CORRECT FOR MHT    
    pList = list(pDict.values()) #uncorrected
    mhtList = fdrcorrection(pList,alpha=0.05)[1] #corrected

    print('\nMHT corrected p-values:')
    for i in range(len(pDict.keys())):
        cell = list(pDict.keys())[i]
        mhtP = mhtList[i]
        print(nameDict[cell],': ',mhtP)        

    
    ##FIGURE 3g - LYMPHOCYTE PD-1 EXPRESSION BY SPATIAL ARCHITECTURE
    
    print('\nFigure 3g') 
       
    print('Calculating PD-1+ counts for each tumor region...')
    pdCounts()
    print("PD-1+ counts calculated. Csvs saved to 'dfCreated' folder.")
     
    #read pd counts file
    df = pd.read_csv(cwd+'/dfCreated/dfPD1PDL1CountsPerc_all.csv', index_col=0)
    
    #read in appropriate mixing score df (either all or avg)
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold ROIs - all ROIs
    coldList = ['pt2post_roi1','pt9post_roi3'] #this is a pre-defined list based on # of immune cells/area
    for roi in coldList:
            dfMix.loc[roi,'spatial'] = 'Cold'
    
    #now add the spatial header into the target column in the dfCount
    df['spatial'] = list(dfMix['spatial'])
    
    pList = [] #for MHT correction later
    
    #separate out mixed and compartmentalized based on spatial column for significance calculation
    dfMixed = df[['Mixed' in row for row in df['spatial']]]
    dfCompart = df[['Compartmentalized' in row for row in df['spatial']]]

    #calculate p value for each cell type and account for MHT
    for cell in ['A_pd1+','B_pd1+','K_pd1+','C_pd1+','D_pd1+','J_pd1+','X_pd1+','E_pd1+','N_pd1+','H_pd1+']:

        #no cells follow normal distr, so use mann whitney U
        pVal = mannwhitneyu(dfMixed.loc[:,cell],dfCompart.loc[:,cell],alternative='greater')[1] #significant p value means mixed is greater than comaprt    
        pList.append(pVal)
    
    #MHT Correction    
    mhtList = fdrcorrection(pList,alpha=0.05)[1]
        
    #create % Lymphocyte PD-1+ figure  
    cellList = ['A_pd1+','B_pd1+','C_pd1+']
    
    dfList = [dfCompart,dfMixed]
        
    #empty list to store two traces for plotting
    traces = []
    
    #loop through the two dfs
    for i in range(len(dfList)):
    
        #get df
        df = dfList[i]
    
        #empty lists to store data for plotting - one per trace
        xAxis = [] #will store label names
        yList = [] #will store count values to plot in boxplot
    
        for col in cellList: #for each column you want to plot
    
            for row in df.index: #for each row of patients
                xAxis.append(nameDict[col[0]]) #add the column's real class name to the xAxis label list (each label should repeat itself nine times)
                yList.append(df.loc[row,col]*100) #add the column/row's value (aka actual value to plot)
    
        #naming of traces
        if i == 0:
            spat = 'Compartmentalized'
        elif i == 1:
            spat = 'Mixed'
    
        #create one trace for pre and one for post
        trace = go.Box(x=xAxis,y=yList,name=spat,boxmean=True,boxpoints='all',marker={'color':colorDict[spat]},line={'color':colorDict[spat]})
        traces.append(trace)
    
    layout = {'xaxis':{'title':'Cell Phenotype'},'yaxis':{'title':'% of Cells PD-1+'},'boxmode':'group'}
    
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure3g.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 3g saved to 'figures' folder.")
    
    print('\nMHT corrected p-values:')
    print('CD8+ T Cell: ',mhtList[0])
    print('CD4+ T Helper: ',mhtList[1])
    print('B Cell: ',mhtList[3])
 
        
    ##FIGURE 3h - APC KI-67 EXPRESSION BY SPATIAL ARCHITECTURE
    
    print('\nFigure 3h') 
    
    #calculate KI-67 positive cell counts
    print('Calculating Ki-67+ cell counts for each tumor region...')
    kiCounts()
    print("Ki-67+ cell counts calculated. Csvs saved to 'dfCreated' folder.")

    #read ki-76 counts csv
    df = pd.read_csv(cwd+'/dfCreated/dfKi67CountsPerc_all.csv', index_col=0)
    
    #read in appropriate mixing score df (either all or avg)
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold ROIs - all ROIs
    coldList = ['pt2post_roi1','pt9post_roi3'] #this is a pre-defined list based on # of immune cells from calculateCold()
    for roi in coldList:
            dfMix.loc[roi,'spatial'] = 'Cold'
    
    #now add the spatial header into the spatial column to the Kisum df
    df['spatial'] = list(dfMix['spatial'])
        
    #separate out mixed and compartmentalized based on spatial column for significance calculation
    dfMixed = df[['Mixed' in row for row in df['spatial']]]
    dfCompart = df[['Compartmentalized' in row for row in df['spatial']]]

    #create figure for APC Ki67+ expression
    traceM = go.Box(y=dfMixed.loc[:,'X_ki67+']*100,name='Mixed',boxpoints='all',marker={'color':'#FF7F03'},line={'color':'#FF7F03'})
    traceC = go.Box(y=dfCompart.loc[:,'X_ki67+']*100,name='Compartmentalized',boxpoints='all',marker={'color':'#2CA02C'},line={'color':'#2CA02C'})
    
    traces = [traceC,traceM]
    layout = {'xaxis':{'title':'Spatial Architecture'},'yaxis':{'title':'% of APCs Ki-67+'},'boxmode':'group'}
    fig = {'data':traces,'layout':layout}
    
    pio.write_image(fig,cwd+'/figures/figure3h.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Figure 3h saved to 'figures' folder.")
    
    #calculate p value for all cell types, then do MHT correction
    pList = [] #empty list to store uncorrected p values
    for cell in ['A_ki67+','B_ki67+','K_ki67+','C_ki67+','D_ki67+','J_ki67+','X_ki67+','E_ki67+','N_ki67+','H_ki67+']:
    
        if cell in ['A_ki67+','J_ki67+','N_ki67+']: #these cells follow normal distribution, so run parametric test
            pTwo = ttest_ind(dfMixed.loc[:,cell],dfCompart.loc[:,cell],equal_var=False)[1] #returns results from TWO-TAIL independent (parametric) t-test        
            pVal = pTwo/2
        
        else: #otherwise run non-parametric mann whitney u for cells that do not follow normal distribution
            pVal = mannwhitneyu(dfMixed.loc[:,cell],dfCompart.loc[:,cell],alternative='less')[1] #significant p value means compart is greater than mixed
    
        pList.append(pVal)
    
    ##CORRECT FOR MHT    
    mhtList = fdrcorrection(pList,alpha=0.05)[1] #corrected
    print('\nMHT corrected p-value:')
    print('Antigen Presenting Cell: ',mhtList[6])
    
        
    ##FIGURE 3i - ASMA MESENCHYMAL CELL DENSITY COLORED BY MIXING SCORE
    
    print('\nFigure 3i')
    
    #read csv with density - this csv gets created in fig1() function
    dfCounts = pd.read_csv(cwd+'/dfCreated/dfDensity_all.csv', index_col=0)
    
    #N = aSMA mesenchymal cells
    count = dfCounts['N']
    
    #sort regions by count in ascending order
    count = count.sort_values()
    
    #get labels of the regions; x axis
    roiLabels = list(count.index)
    
    #get counts to plot; y axis
    countValues = list(count.values)
    
    #read csv with mixing scores
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv',index_col=0)
    
    #get just the mix score column
    mix = dfMix['mixScore']
    
    #empty list to store colors according to mixing score
    colorList = []
    mixLabel = []
    
    #loop through each region and assign color based on its mixing score and immune count (for cold only)
    for roi in roiLabels:
    
        #if not cold
        if roi != 'pt2post_roi1' and roi != 'pt9post_roi3': #determined by calculateCold() function
    
            #compartmentalized (median of averaged primary mixing scores)
            if mix[roi] < 0.107:
                colorList.append('#2CA02C')
                mixLabel.append('Compartmentalized')
            #mixed
            else:
                colorList.append('#FF7F03')
                mixLabel.append('Mixed')
        #cold
        else:
            colorList.append('#E377C2')
            mixLabel.append('Cold')
    
    #make figure
    trace = [go.Bar(x=roiLabels,y=countValues,marker={'color':colorList})]
    
    #make annotations for mixing labels to appear on plot
    annotC = dict(x=37.6,y=500,text='Compartmentalized',showarrow=False,font={'color':'#2CA02C'})
    annotM = dict(x=35.3,y=450,text='Mixed',showarrow=False,font={'color':'#FF7F03'})
    annotC2 = dict(x=35,y=400,text='Cold',showarrow=False,font={'color':'#E377C2'})
    annotList = [annotC,annotC2,annotM]
    
    layout = {'xaxis':{'title':'Tumor Region','automargin':True},'yaxis':{'title':'⍺SMA+ Mesenchymal Cell Density'},'annotations':annotList}
    fig = {'data':trace,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure3i.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 3i saved to 'figures' folder.")
    


def fig4():
    
    '''
    αSMA+ mesenchymal cellular neighborhood clustering.
    This function creates figures 4b-4d, 4f-4g. 
    All subplots are saved to the 'figures' folder.
        
    Input parameters:
        None
        
    Outputs:
        Saves csvs to 'dfCreated' folder for files containing cellular neighborhood data used to generate figures.
        Saves plots to 'figures' folder for figures 4b-4d, 4f-4g 
    '''
  
    #import packages    
    import pandas as pd
    import os
    import plotly.io as pio
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    import plotly.graph_objs as go
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    #get dictionary to rename cells - used for subfigures
    nameDict = cellClassRename()
    
    #get current directory - used to save subfigures
    cwd = os.getcwd()
    
    #predefined color scheme; one distinct color per 11 cell types
    cellColors = {'A':'rgb(127,60,141)','B':'rgb(17,165,121)','K':'rgb(57,105,172)','C':'rgb(242,183,1)','D':'rgb(231,63,116)','J':'rgb(128,186,90)','X':'rgb(230,131,16)','E':'rgb(0,134,149)','N':'rgb(207,28,144)','F':'rgb(249,123,114)','H':'rgb(165,170,153)'}

    #change numbering of columns from 0-6 to 1-7
    colDict = {'clust0PercList':'1','clust1PercList':'2','clust2PercList':'3','clust3PercList':'4','clust4PercList':'5','clust5PercList':'6','clust6PercList':'7'}


    ##FIGURE 4b - AVERAGE αSMA+ MESENCHYMAL CELL NEIGHBORHOOD CLUSTER COMPOSITION

    print('\nFigure 4b')

    #make neighborhoods for αSMA+ seed cells (N), distThresh = 30µm = 60px
    print('Creating αSMA+ cellular neighborhoods...')
    makeNeighborhoods('N',60)
    print("Neighborhoods created. Csv saved to 'dfCreated' folder.")

    #optionally run elbowMethod() function to determine optimal number of clusters
    print('Running Elbow Method to determine optimal number of clusters...')
    elbowMethod('dfNeighborhoodClusterN60',15) #move forward with k=7 based on these results
    
    #cluster αSMA+ cellular neighborhoods based on composition, k=7
    print('Clustering αSMA+ cellular neighborhoods...')
    clusterNeighborhoods('dfNeighborhoodClusterN60',7)
    print("Neighborhoods clustered. Csv saved to 'dfCreated' folder.")
        
    #read csv
    df = pd.read_csv(cwd+'/dfCreated/dfNeighClusteredN60k7.csv', index_col=0)
    
    #filter df to not include the index or file columns
    df = df.drop(['index','file'],axis=1)
    
    #groupby cluster column and take the averages of all of the other columns for each group
    dfCluster = df.groupby(['cluster']).mean()
    
    #transpose df so columns are clusters and rows are cell types
    dfClusterT = dfCluster.T
    
    #rename clusters to be 1-6 rather than 0-5
    dfClusterT = dfClusterT.rename(columns={0:'1',1:'2',2:'3',3:'4',4:'5',5:'6',6:'7'})
    
    #re-order the rows of dfClusterT such that immune cells are first, followed by stromal, followed by tumor
    dfClusterT = dfClusterT.reindex(["countA%", "countB%", "countK%",'countC%','countD%','countJ%','countX%','countE%','countN%','countH%'])
    
    #create x axis labels
    xAxis = list(dfClusterT.columns)    
    
    #create one trace per cluster (aka per column in dfClusterT)
    traces = []
    for i in range(len(dfClusterT)):
        #get cell name from A-X label using nameDict
        cell = nameDict[dfClusterT.index[i][-2]]
    
        #create trace and add to traces
        trace = {'type':'bar','x':xAxis,'y':list(dfClusterT.iloc[i,0:]),'name':cell,'marker':{'color':cellColors[dfClusterT.index[i][-2]]}}
        traces.append(trace)
    
    #plot all clusters along x, and have one bar with all cell types (STACKED)
    layout = {'xaxis':{'title':'Cluster'},'yaxis':{'title':'Fraction Present'},'barmode':'stack'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure4b.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 4b saved to 'figures' folder.")


    ##FIGURE 4c - HIERARCHICAL CLUSTERING OF ASMA NEIGHBORHOOD CLUSTERS

    print('\nFigure 4C')
    
    #calculate abundance of each cluster in each tumor region, save data to csvs
    print('Calculating counts of clusters per tumor region...')
    clusterCountPerROI('dfNeighClusteredN60k7')
    print("Counts calculated. Csvs saved to 'dfCreated' folder.")

    #load df
    df = pd.read_csv(cwd+'/dfCreated/dfClustCountsN60k7_avg.csv',index_col=0)
    
    #filter to just primaries
    df = df[['pre' in pt for pt in df.index]]
        
    #add a new column for one color - patient number (1-9)
    df['patient'] = df.index.str.slice(2,3)
    
    #make color vector dictionary; two in total for two color palettes
    palette1 = dict(zip(df.patient.unique(),sns.color_palette('tab10',n_colors = 9))) ##9 colors for 9 patients
    #make groupings
    grouping1 = df.patient #group off of pt column
    #map palettes onto the corresponding grouping series; names will be used as labels for color columns
    ptColors1 = pd.Series(grouping1,name='Patient ID').map(palette1)
    dfColors = pd.concat([ptColors1], axis=1)
    
    #get only the perc columns
    df = df[[col for col in df.columns if 'Perc' in col]]
    df = df.rename(mapper=colDict, axis=1)
    
    #normalize to log10+1 scale of % out of 100
    dfLog = np.log10((df*100)+1)
    
    #plot the clustermap with the colors
    graph = sns.clustermap(dfLog,method='ward',metric='euclidean',cmap='vlag',row_colors=dfColors,yticklabels=True) #yticklabels = true to show all patient labels; row_colors should be on the values of ptColors
    ax = graph.ax_heatmap
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Patient")
    #manually add colors for the two "groups" in illustrator
    bottom,top = ax.get_ylim()
    ax.set_ylim(bottom+.5,top-.5)

    plt.savefig(cwd+'/figures/figure4c.png',format='png',bbox_inches='tight') #saves figure to 'figures' folder as a png file
    plt.close()
    
    print("Figure 4c saved to 'figures' folder.")


    ##FIGURE 4d - AVERAGE ASMA NEIGHBORHOOD CLUSTER COMPOSITION OF TWO GROUPS

    print('\nFigure 4d')

    #read csv
    df = pd.read_csv(cwd+'/dfCreated/dfClustCountsN60k7_avg.csv', index_col=0)
    
    #filter to primary only
    df = df[['pre' in pt for pt in df.index]]
       
    #filter to just columns with Perc in name
    df = df[df.columns[['Perc' in col for col in list(df.columns)]]]
    
    #add group assignment
    group1 = ['pt1pre','pt3pre','pt6pre','pt8pre'] #based off clustering analysis
    df['group'] = np.where(df.index.isin(group1), '1','2')
    
    #groupby cluster column and take the averages of all of the other columns for each group
    dfGroup = df.groupby(['group']).mean()
    dfGroup
    
    xAxis = list(dfGroup.index)
    traces = []
    
    for col in list(dfGroup.columns):
        trace = go.Bar(x=xAxis,y=dfGroup[col],name=colDict[col])
        traces.append(trace)
    
    layout = {'barmode':'stack','xaxis':{'title':'Group'},'yaxis':{'title':'Fraction Present'}}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure4d.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 4d saved to 'figures' folder.")


    ##FIGURE 4f - ASMA NEIGHBORHOOD GROUP SURVIVAL ANALYSIS

    print('\nFigure 4f')

    #get clinical df data
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv', index_col=0)
    dfClin.set_index('sample',inplace=True) #set index labels to be the patient column
    
    #filter df so only desired tumors are analyzed
    dfClin = dfClin[['pre' in row for row in dfClin.index]]
    dfClin.index = dfClin.index.str.slice(2,3)
    dfClin['dtr'] = pd.to_numeric(dfClin['dtr'],errors='coerce')
    
    dfClinAvg = dfClin.groupby(['sample']).mean()
    
    #pre-defined groups based off fig 4c
    dfLow = dfClinAvg[dfClinAvg.index.isin(['1','3','6','8'])] #group1
    dfHigh = dfClinAvg[dfClinAvg.index.isin(['2','4','5','7','9'])] #group2
            
    #fit KM model
    kmf = KaplanMeierFitter()
    
    #create subplot to hold both curves
    plt.figure() #new figure for each call to this function
    ax = plt.subplot(111)
    
    #dfLow curve
    kmf.fit(dfLow['dtr'], event_observed=[1]*len(dfLow), label="Group 1")
    kmf.plot(ax=ax,ci_show=False) #can change to True to see confidence intervals
    
    #dfHigh curve
    kmf.fit(dfHigh['dtr'], event_observed=[1]*len(dfHigh), label="Group 2")
    kmf.plot(ax=ax,ci_show=False)
    
    plt.ylim(0, 1) #set y axis to go from zero to one
    plt.margins(0) #make (0,0) start in bottom left corner rather than a little bit off
    
    #run logrank test for p-value calculation
    results = logrank_test(dfLow['dtr'],dfHigh['dtr'],[1]*len(dfLow),[1]*len(dfHigh),alpha=0.95)
    
    #title and axis labels
    plt.ylabel('Progression Free Survival (%)')
    plt.xlabel('Days')
    plt.text(x=1150,y=0.7,s='p = '+str(round(results.p_value,3)))
    
    plt.savefig(cwd+'/figures/figure4f.png',format='png') #saves figure to 'figures' folder as a png file
    plt.close()
    print("Figure 4f saved to 'figures' folder.")        


    ##FIGURE 4g - PRIMARY TUMOR ASMA NEIGHBORHOOD CLUSTER COMPOSITION

    print('\nFigure 4g')

    #get df from saved csv
    df = pd.read_csv(cwd+'/dfCreated/dfClustCountsN60k7_avg.csv', index_col=0)
    df = df[['pre' in pt for pt in df.index]]
       
    #generate column list to cluster on based on if there is a % in the column name
    dfPerc = df[df.columns[['Perc' in col for col in list(df.columns)]]]
    dfPerc.index = dfPerc.index.str.slice(2,3)
    
    #sort descending cluster 1
    dfPercSort = dfPerc.sort_values('clust0PercList',ascending=False)
    
    xAxis = list(dfPercSort.index)
    traces = []
    
    for col in list(dfPercSort.columns):
        trace = go.Bar(x=xAxis,y=dfPercSort[col],name=colDict[col])
        traces.append(trace)
    
    layout = {'barmode':'stack','xaxis':{'title':'Primary Tumor'},'yaxis':{'title':'Fraction Present'}}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure4g.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 4g saved to 'figures' folder.")
    
    
    
def table3():
    
    '''
    This function creates table 3, the coefficient of variation across cell types in the cohort.
    
    Note: This function assumes densityToCsv() has already been run (this happens in fig1() function).
    
    Input parameters:
        None
    
    Outputs:
        Saves table to 'figures' folder for table 3
        
    '''
    
    ##TABLE 3 - AVERAGE COEFFICIENT OF VARIATION PER CELL TYPE
    #import packages    
    import pandas as pd
    import os

    print('\nTable 3')
    
    nameDict = cellClassRename()
    
    cwd = os.getcwd()
 
    #get df of densities
    df = pd.read_csv(cwd+'/dfCreated/dfDensity_all.csv',index_col=0)
    
    #empty dict to store stats values across all tumors
    tumCvDict = {}
    
    #subset df into smaller dfs per 18 tumors
    for tum in df.index.str.slice(0,6).unique():
        
        dfTumor = df[df.index.str.slice(0,6) == tum].iloc[:,1:] #add iloc to drop dtr column
    
        #calculate range for each cell type
        rng = dfTumor.max() - dfTumor.min()
        
        if rng.all() != 0: #drop the 2 tumors (pt2pos,pt6pos) with only one region
            tumCvDict[tum] = (dfTumor.std())/(dfTumor.mean()) #coefficient of variation = std/mean
    
    cvAvgSeries = round(pd.DataFrame(tumCvDict).T.mean(),3)
    dfCV = pd.DataFrame(cvAvgSeries,columns=['Average Coefficient of Variation']).sort_values('Average Coefficient of Variation',ascending=False)
    dfCV.index = dfCV.index.str.slice(0,1)
    dfCV = dfCV.rename(index=nameDict)
    dfCV.to_csv(cwd+'/figures/table3.csv')

    print("Table 3 saved to 'figures' folder.")

    

def supp1():
    
    '''
    This function creates supplemental figure 1c.
    Subplot is saved to the 'figures' folder.
    
    Note: This function assumes densityToCsv() has already been run (this happens in fig1() function), as it relies on density csvs to generate figure 1c.
    
    Input parameters:
        None
        
    Outputs:
        Saves plot to 'figures' folder for figure 1c
    '''
  
    import pandas as pd
    import numpy as np
    from scipy.stats import entropy,f_oneway
    import plotly.graph_objs as go
    import plotly.io as pio
    import os
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
  
    ##SUPPLEMENTAL FIGURE 1c - KL DIVERGENCE, IMMUNE ONLY
    
    print('\nSupplemental Figure 1c')
    
    #get current directory
    cwd = os.getcwd()
    
    #rename xAxis variables to 1P, 1R, etc
    xDictRename = {'pre':'P','pos':'R'}

    #get df of cell counts - using the perc composition with ALL cells (immune, tumor, aSMA)
    dfAll = pd.read_csv(cwd+'/dfCreated/dfDensity_all.csv',index_col=0)
    dfAvg = pd.read_csv(cwd+'/dfCreated/dfDensity_avg.csv',index_col=0)
    
    #get average of dfAvg - average values per patient (P+R)/2
    dfAvgCopy = dfAvg.copy()
    dfAvgCopy.index = dfAvgCopy.index.str.slice(2,3)
    dfAvgAvg = dfAvgCopy.groupby(['patient']).mean()
    
    #get average of primary only and recurrent only 
    dfPrim = dfAll[['pre' in i for i in dfAll.index]]
    dfPrimAvg = dfPrim.mean()
    
    dfRec = dfAll[['pos' in i for i in dfAll.index]]
    dfRecAvg = dfRec.mean()
    
    #get average of cohort
    dfCohortAvg = dfAvgAvg.mean()
    
    #get clinical df with anatomical location to add a column for site on dfAll
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv',index_col=0)
    dfClin.set_index('sample',inplace=True)
    
    dfAllCopy = dfAll.copy()
    dfAllCopy['anatomy'] = dfClin['anatomy']
    anatDict = dict(zip(dfAllCopy.index, dfAllCopy.anatomy)) #store roi:anatomy pairs for later use
    
    #get average of each anatomy region (OC,OP,L)
    dfOCAvg = dfAllCopy[dfAllCopy['anatomy'] == 'oral cavity'].mean()
    dfOPAvg = dfAllCopy[dfAllCopy['anatomy'] == 'oropharynx'].mean()
    dfLAvg = dfAllCopy[dfAllCopy['anatomy'] == 'larynx'].mean()
    
    #empty list to store KL divergence
    entListA = [] #for intra-tumoral, compares to single tumor's average
    entListB = [] #for intra-patient, compares to single patient's average (P+R)/2
    entListC = [] #for inter-patient, compares to all primary, or all recurrent depending on which one the region is from
    entListD = [] #for inter-patient, compares to cohort average across ALL pts and tumors, regardless of P or R
    entListY = [] #for intra-site, compares to same anatomic site's average across the cohort

    #create list of rois (patient samples)
    ptList = list(dfAll.index)
    
    #loop through ROIs
    for pt in ptList:
        #get row values for A-E counts for each ROI (per the pt variable); immune ONLY!
        countRoi = dfAll.loc[pt,'A':'E'] 
        
        #get character(s) to match other df's index column
        patA = pt[0:6] #same tumor
        patB = pt[2] #same pt
        
        #get row values for A-H counts for each tumor (averaged, use the pat variables)
        countAvgA = dfAvg.loc[patA,'A':'E']
        countAvgB = dfAvgAvg.loc[patB,'A':'E']
        
        if 'pre' in pt:
            countAvgC = dfPrimAvg['A':'E']
            
        elif 'pos' in pt:
            countAvgC = dfRecAvg['A':'E']
    
        #values for matching anatomic site
        if anatDict[pt] == 'oral cavity':
            countAvgY = dfOCAvg['A':'E']
        elif anatDict[pt] == 'oropharynx':
            countAvgY = dfOPAvg['A':'E']
        elif anatDict[pt] == 'larynx':
            countAvgY = dfLAvg['A':'E']
    
        #cohort avg
        countAvgD = dfCohortAvg['A':'E']
        
        #calculate the KL Divergence using the entropy function (log2) for each avg distribution; sample distribution's divergence from the tumor's, patient's, cohort's average (qk) distribution
        entA = entropy(countRoi,qk=countAvgA,base=2) #intra-tumor
        entB = entropy(countRoi,qk=countAvgB,base=2) #intra-patient
        entC = entropy(countRoi,qk=countAvgC,base=2) #inter-patient (p or r)
        entD = entropy(countRoi,qk=countAvgD,base=2) #inter-patient (all)
        entY = entropy(countRoi,qk=countAvgY,base=2) #intra-anatomic site

        #store KL divergence in list
        entListA.append(entA)
        entListB.append(entB)
        entListC.append(entC)
        entListD.append(entD)
        entListY.append(entY)
    
    #create new df with KL divergence info stored in one column
    dfEnt = pd.DataFrame({'Intra-Tumor (P or R only)':entListA,'Intra-Patient (P and R)':entListB,'Inter-Patient (P or R only)':entListC,'Inter-Patient (Same Anatomic Site)':entListY,'Inter-Patient (all)':entListD},index=ptList)
    
    #empty list to store x axis labels
    xAxis = []
    #rename xAxis variables to 1P, 1R, etc
    for x in dfEnt.index:
        xLabel = x[2]+xDictRename[x[3:6]]
        xAxis.append(xLabel) #add label to list
    
    traces = []
    for col in list(dfEnt.columns):        
        trace = go.Box(x=xAxis,y=dfEnt[col],boxpoints='all',name=col)
        traces.append(trace)
    layout = {'xaxis':{'title':'Tumor'},'yaxis':{'title':'KL Divergence'},'boxmode':'group'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/supp1c.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Supplemental Figure 1c saved to 'figures' folder.")

    #run ANOVA to see if means are the same or different; p is significant, so not all group means are equal! 
    anP = f_oneway(entListA,entListB,entListC,entListY,entListD)[1]
    print('One way ANOVA p-value: ',anP) #one way anova

    #if ANOVA is significant, run Tukey post hoc analysis to determine which pairs are significant
    if anP < 0.05:        
        dfTukey = pd.DataFrame({'value':entListA+entListB+entListC+entListY+entListD,'label':np.repeat(list(dfEnt.columns),repeats=len(dfEnt))})
        tukey = pairwise_tukeyhsd(endog=dfTukey['value'],groups=dfTukey['label'],alpha=0.05)
        print('Tukey HSD Post-Hoc Test: ',tukey)



def supp2():
    
    '''
    This function creates supplemental figures 2a-2e.
    All subplots are saved to the 'figures' folder.
    
    Note: This function assumes calcluateMixing(), pdCounts(), and kiCounts() have already been run (this happens in fig3() function).
    
    Input parameters:
        None
        
    Outputs:
        Saves plots to 'figures' folder for figures 2a-2e
    '''
  
    import pandas as pd
    import numpy as np
    from scipy.stats import f_oneway
    import plotly.graph_objs as go
    import plotly.io as pio
    import os
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
  
    #get current directory
    cwd = os.getcwd()
    
    #rename xAxis variables to 1P, 1R, etc
    xDictRename = {'pre':'P','pos':'R'}

    #get mixed/compart designation
    def map_mix(mix):
    
        if mix < 0.107:
            return 'Compartmentalized'
        elif mix > 0.107:
            return 'Mixed'


    ##SUPPLEMENTAL FIGURE 2a - ABSOLUTE DIFFERENCE IN MIXING SCORE
    
    print('\nSupplemental Figure 2a')

    #get df of mixing scores
    dfAll = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv',index_col=0)
    dfAvg = pd.read_csv(cwd+'/dfCreated/dfMixingScores_avg.csv',index_col=0)
    
    #get average of dfAvg - average values per patient (P+R)/2
    dfAvgCopy = dfAvg.copy()
    dfAvgCopy.index = dfAvgCopy.index.str.slice(2,3)
    dfAvgAvg = dfAvgCopy.groupby(['patient']).mean()
    
    #get average of primary only and recurrent only 
    dfPrim = dfAll[['pre' in i for i in dfAll.index]]
    dfPrimAvg = dfPrim.mean()
    
    dfRec = dfAll[['pos' in i for i in dfAll.index]]
    dfRecAvg = dfRec.mean()
    
    #get average of cohort
    dfCohortAvg = dfAvgAvg.mean()
    
    #get clinical df with anatomical location to add a column for site on dfAll
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv',index_col=0)
    dfClin.set_index('sample',inplace=True)
    
    dfAllCopy = dfAll.copy()
    dfAllCopy['anatomy'] = dfClin['anatomy']
    anatDict = dict(zip(dfAllCopy.index, dfAllCopy.anatomy)) #store roi:anatomy pairs for later use
    
    #get average of each anatomy region (OC,OP,L)
    dfOCAvg = dfAllCopy[dfAllCopy['anatomy'] == 'oral cavity'].mean()
    dfOPAvg = dfAllCopy[dfAllCopy['anatomy'] == 'oropharynx'].mean()
    dfLAvg = dfAllCopy[dfAllCopy['anatomy'] == 'larynx'].mean()
    
    #empty list to store differences
    entListA = [] #for intra-tumoral, compares to single tumor's average
    entListB = [] #for intra-patient, compares to single patient's average (P+R)/2
    entListC = [] #for inter-patient, compares to all primary, or all recurrent depending on which one the region is from
    entListD = [] #for inter-patient, compares to cohort average across ALL pts and tumors, regardless of P or R
    entListY = [] #for intra-site, comapres to same anatomic site's average across the cohort
    
    #create list of rois (patient samples)
    ptList = list(dfAll.index)
    
    #loop through ROIs
    for pt in ptList:
        #get mixing scores for each ROI (per the pt variable)
        countRoi = dfAll.loc[pt,'mixScore'] 
        
        #get character(s) to match other df's index column
        patA = pt[0:6]
        patB = pt[2]
        
        #get mixing scores for each tumor (averaged, use the pat variables)
        countAvgA = dfAvg.loc[patA,'mixScore']
        countAvgB = dfAvgAvg.loc[patB,'mixScore']
        
        if 'pre' in pt:
            countAvgC = dfPrimAvg['mixScore']
            
        elif 'pos' in pt:
            countAvgC = dfRecAvg['mixScore']
        
        #values for matching anatomic site
        if anatDict[pt] == 'oral cavity':
            countAvgY = dfOCAvg['mixScore']
        elif anatDict[pt] == 'oropharynx':
            countAvgY = dfOPAvg['mixScore']
        elif anatDict[pt] == 'larynx':
            countAvgY = dfLAvg['mixScore']
    
        countAvgD = dfCohortAvg['mixScore']
        
        #calculate the difference between the region's mixing score and the average mixing score at 5 levels the tumor's, patient's, cohort's average (qk) distribution
        entA = abs(countRoi-countAvgA) #intra-tumor
        entB = abs(countRoi-countAvgB) #intra-patient
        entC = abs(countRoi-countAvgC) #inter-patient (p or r)
        entD = abs(countRoi-countAvgD) #inter-patient (all)
        entY = abs(countRoi-countAvgY) #intra-anatomic site
    
        #store diffs in list
        entListA.append(entA)
        entListB.append(entB)
        entListC.append(entC)
        entListD.append(entD)
        entListY.append(entY)
        
    #create new df with diff info stored in one column
    dfEnt = pd.DataFrame({'Intra-Tumor (P or R only)':entListA,'Intra-Patient (P and R)':entListB,'Inter-Patient (P or R only)':entListC,'Inter-Patient (Same Anatomic Site)':entListY,'Inter-Patient (all)':entListD},index=ptList)
    
    #empty list to store x axis labels
    xAxis = []
    #rename xAxis variables to 1P, 1R, etc
    xDictRename = {'pre':'P','pos':'R'}
    for x in dfEnt.index:
        xLabel = x[2]+xDictRename[x[3:6]]
        xAxis.append(xLabel) #add label to list
    
    traces = []
    for col in list(dfEnt.columns):
        trace = go.Box(x=xAxis,y=dfEnt[col],boxpoints='all',name=col)
        traces.append(trace)
    
    layout = {'xaxis':{'title':'Tumor'},'yaxis':{'title':'Absolute Difference'},'boxmode':'group'}
    fig = {'data':traces,'layout':layout}
    
    pio.write_image(fig,cwd+'/figures/supp2a.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    print("Supplemental Figure 2a saved to 'figures' folder.")

    #run ANOVA to see if means are the same or different; p is significant, so not all group means are equal! 
    anP = f_oneway(entListA,entListB,entListC,entListY,entListD)[1]
    print('One way ANOVA p-value: ',anP) #one way anova

    #if ANOVA is significant, run Tukey post hoc analysis to determine which pairs are significant
    if anP < 0.05:        
        dfTukey = pd.DataFrame({'value':entListA+entListB+entListC+entListY+entListD,'label':np.repeat(list(dfEnt.columns),repeats=len(dfEnt))})
        tukey = pairwise_tukeyhsd(endog=dfTukey['value'],groups=dfTukey['label'],alpha=0.05)
        print('Tukey HSD Post-Hoc Test: ',tukey)


    ##SUPPLEMENTAL FIGURE 2b - PD-1+ CD8+ T CELL BOOTSTRAPPING
    
    print('\nSupplemental Figure 2b')

    #get pd1 counts df
    df = pd.read_csv(cwd+'/dfCreated/dfPD1PDL1CountsPerc_all.csv', index_col=0)
    
    #filter to just cd8s
    df = df[['A_pd1+']]
    
    #read in mix df
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold ROIs - all ROIs
    coldList = ['pt2post_roi1','pt9post_roi3'] #this is a pre-defined list based on # of immune cells
    for roi in coldList:
            dfMix.loc[roi,'spatial'] = 'Cold'
    
    #now add the spatial header into the spatial column to the df
    df['spatial'] = list(dfMix['spatial'])
    
    #empty df to store the bootstrapped data
    dfRandom = pd.DataFrame()
    
    #empty lists to store means of each group
    mixMeanList = []
    compMeanList = []
    
    #bootstrap through 100 times
    for i in range(100):
        
        for tum in df.index.str.slice(0,6).unique(): #for each of the 18 tumors    
            #subset df to only include unique tumors
            dfTumor = df[df.index.str.slice(0,6) == tum]
    
            #choose a random roi for that tumor
            row = np.random.choice(dfTumor.index,size=1)
    
            #add the randomly selected row to dfRandom
            dfRandom = dfRandom.append(dfTumor.loc[row,:])
    
        #separate out mixed and compartmentalized based on spatial column
        dfMixed = dfRandom[['Mixed' in row for row in dfRandom['spatial']]]
        dfCompart = dfRandom[['Compartmentalized' in row for row in dfRandom['spatial']]]
    
        #calculate mean %pd-1+ values for each cell for mixed and compart groups and store in lists - get one mean per 100 bootstraps through each of the 18 tumors
        mixMeanList.append(dfMixed.mean()[0]*100)
        compMeanList.append(dfCompart.mean()[0]*100)

    traceC = go.Histogram(x=compMeanList,name='Compartmentalized',nbinsx=100,marker_color='#2CA02C')
    traceM = go.Histogram(x=mixMeanList,name='Mixed',nbinsx=100,marker_color='#FF7F03')
    traces = [traceC,traceM]
    
    layout = {'xaxis':{'title':'% PD-1+ CD8+ T Cells'},'yaxis':{'title':'Frequency'}}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/supp2b.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Supplemental Figure 2b saved to 'figures' folder.")


    ##SUPPLEMENTAL FIGURE 2c - PD-1+ CD4+ T HELPER CELL BOOTSTRAPPING
    
    print('\nSupplemental Figure 2c')

    #get pd1 counts df
    df = pd.read_csv(cwd+'/dfCreated/dfPD1PDL1CountsPerc_all.csv', index_col=0)
    
    #filter to just t helpers
    df = df[['B_pd1+']]
    
    #read in mix df
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold ROIs - all ROIs
    coldList = ['pt2post_roi1','pt9post_roi3'] #this is a pre-defined list based on # of immune cells
    for roi in coldList:
            dfMix.loc[roi,'spatial'] = 'Cold'
    
    #now add the spatial header into the spatial column to the df
    df['spatial'] = list(dfMix['spatial'])
    
    #empty df to store the bootstrapped data
    dfRandom = pd.DataFrame()
    
    #empty lists to store means of each group
    mixMeanList = []
    compMeanList = []
    
    #bootstrap through 100 times
    for i in range(100):
        
        for tum in df.index.str.slice(0,6).unique(): #for each of the 18 tumors
    
            #subset dfPd to only include unique tumors
            dfTumor = df[df.index.str.slice(0,6) == tum]
    
            #choose a random roi for that tumor
            row = np.random.choice(dfTumor.index,size=1)
    
            #add the randomly selected row to dfRandom
            dfRandom = dfRandom.append(dfTumor.loc[row,:])
    
        #separate out mixed and compartmentalized based on spatial column
        dfMixed = dfRandom[['Mixed' in row for row in dfRandom['spatial']]]
        dfCompart = dfRandom[['Compartmentalized' in row for row in dfRandom['spatial']]]
    
        #calculate mean %pd-1+ values for each cell for mixed and compart groups and store in lists - get one mean per 100 bootstraps through each of the 18 tumors
        mixMeanList.append(dfMixed.mean()[0]*100)
        compMeanList.append(dfCompart.mean()[0]*100)

    traceC = go.Histogram(x=compMeanList,name='Compartmentalized',nbinsx=100,marker_color='#2CA02C')
    traceM = go.Histogram(x=mixMeanList,name='Mixed',nbinsx=100,marker_color='#FF7F03')
    traces = [traceC,traceM]
    
    layout = {'xaxis':{'title':'% PD-1+ CD4+ T Helper Cells'},'yaxis':{'title':'Frequency'}}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/supp2c.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Supplemental Figure 2c saved to 'figures' folder.")


    ##SUPPLEMENTAL FIGURE 2d - PD-1+ B CELL BOOTSTRAPPING
    
    print('\nSupplemental Figure 2d')

    #get pd1 counts df
    df = pd.read_csv(cwd+'/dfCreated/dfPD1PDL1CountsPerc_all.csv', index_col=0)
    
    #filter to just b cells
    df = df[['C_pd1+']]
    
    #read in mix df
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold ROIs - all ROIs
    coldList = ['pt2post_roi1','pt9post_roi3'] #this is a pre-defined list based on # of immune cells
    for roi in coldList:
            dfMix.loc[roi,'spatial'] = 'Cold'
    
    #now add the spatial header into the spatial column to the df
    df['spatial'] = list(dfMix['spatial'])
    
    #empty df to store the bootstrapped data
    dfRandom = pd.DataFrame()
    
    #empty lists to store means of each group
    mixMeanList = []
    compMeanList = []
    
    #bootstrap through 100 times
    for i in range(100):
        
        for tum in df.index.str.slice(0,6).unique(): #for each of the 18 tumors
    
            #subset dfPd to only include unique tumors
            dfTumor = df[df.index.str.slice(0,6) == tum]
    
            #choose a random roi for that tumor
            row = np.random.choice(dfTumor.index,size=1)
    
            #add the randomly selected row to dfRandom
            dfRandom = dfRandom.append(dfTumor.loc[row,:])
    
        #separate out mixed and compartmentalized based on spatial column
        dfMixed = dfRandom[['Mixed' in row for row in dfRandom['spatial']]]
        dfCompart = dfRandom[['Compartmentalized' in row for row in dfRandom['spatial']]]
    
        #calculate mean %pd-1+ values for each cell for mixed and compart groups and store in lists - get one mean per 100 bootstraps through each of the 18 tumors
        mixMeanList.append(dfMixed.mean()[0]*100)
        compMeanList.append(dfCompart.mean()[0]*100)

    traceC = go.Histogram(x=compMeanList,name='Compartmentalized',nbinsx=100,marker_color='#2CA02C')
    traceM = go.Histogram(x=mixMeanList,name='Mixed',nbinsx=100,marker_color='#FF7F03')
    traces = [traceC,traceM]
    
    layout = {'xaxis':{'title':'% PD-1+ B Cells'},'yaxis':{'title':'Frequency'}}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/supp2d.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    print("Supplemental Figure 2d saved to 'figures' folder.")


    ##SUPPLEMENTAL FIGURE 2e - Ki-67+ APC BOOTSTRAPPING
    
    print('\nSupplemental Figure 2e')

    #get ki67 counts df
    df = pd.read_csv(cwd+'/dfCreated/dfKI67CountsPerc_all.csv', index_col=0)
    
    #filter to just apcs
    df = df[['X_ki67+']]
    
    #read in mix df
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold ROIs - all ROIs
    coldList = ['pt2post_roi1','pt9post_roi3'] #this is a pre-defined list based on # of immune cells
    for roi in coldList:
            dfMix.loc[roi,'spatial'] = 'Cold'
    
    #now add the spatial header into the spatial column to the df
    df['spatial'] = list(dfMix['spatial'])
    
    #empty df to store the bootstrapped data
    dfRandom = pd.DataFrame()
    
    #empty lists to store means of each group
    mixMeanList = []
    compMeanList = []
    
    #bootstrap through 100 times
    for i in range(100):
        
        for tum in df.index.str.slice(0,6).unique(): #for each of the 18 tumors
    
            #subset dfPd to only include unique tumors
            dfTumor = df[df.index.str.slice(0,6) == tum]
    
            #choose a random roi for that tumor
            row = np.random.choice(dfTumor.index,size=1)
    
            #add the randomly selected row to dfRandom
            dfRandom = dfRandom.append(dfTumor.loc[row,:])
    
        #separate out mixed and compartmentalized based on spatial column
        dfMixed = dfRandom[['Mixed' in row for row in dfRandom['spatial']]]
        dfCompart = dfRandom[['Compartmentalized' in row for row in dfRandom['spatial']]]
    
        #calculate mean %ki-67+ values for each cell for mixed and compart groups and store in lists - get one mean per 100 bootstraps through each of the 18 tumors
        mixMeanList.append(dfMixed.mean()[0]*100)
        compMeanList.append(dfCompart.mean()[0]*100)

    traceC = go.Histogram(x=compMeanList,name='Compartmentalized',nbinsx=100,marker_color='#2CA02C')
    traceM = go.Histogram(x=mixMeanList,name='Mixed',nbinsx=100,marker_color='#FF7F03')
    traces = [traceC,traceM]
    
    layout = {'xaxis':{'title':'% Ki-67+ APCs'},'yaxis':{'title':'Frequency'}}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/supp2e.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    print("Supplemental Figure 2e saved to 'figures' folder.")



def supp3():
    
    '''
    This function creates supplemental figures 3b-3f.
    All subplots are saved to the 'figures' folder.
    
    Note: This function assumes all neighborhood clustering functions have been run (this happens in fig4() function).
    
    Input parameters:
        None
        
    Outputs:
        Saves plots to 'figures' folder for figures 3b-3f
    '''
  
    import pandas as pd
    import plotly.graph_objs as go
    import plotly.io as pio
    import os  
    
    #get current directory
    cwd = os.getcwd()
    
    #rename xAxis variables to 1P, 1R, etc
    xDictRename = {'pre':'P','pos':'R'}

    #get mixed/compart designation
    def map_mix(mix):
    
        if mix < 0.107:
            return 'Compartmentalized'
        elif mix > 0.107:
            return 'Mixed'

    colDict = {'clust0PercList':'1','clust1PercList':'2','clust2PercList':'3','clust3PercList':'4','clust4PercList':'5','clust5PercList':'6','clust6PercList':'7'}


    ##SUPPLEMENTAL FIGURE 3b - ASMA CLUSTER PRESENCE ACROSS TUMOR REGIONS
    
    print('\nSupplemental Figure 3b')

    #read csv
    df = pd.read_csv(cwd+'/dfCreated/dfClustCountsN60k7_all.csv', index_col=0)
    df = df[df.columns[['Perc' in col for col in list(df.columns)]]]
    
    xAxis = list(df.index)
    traces = []
    
    for col in list(df.columns):
        trace = go.Bar(x=xAxis,y=df[col],name=colDict[col])
        traces.append(trace)
    
    layout = {'xaxis':{'title':'Tumor Region'},'yaxis':{'title':'Fraction Present'},'barmode':'stack'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/supp3b.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Supplemental Figure 3b saved to 'figures' folder.")


    ##SUPPLEMENTAL FIGURE 3c - ASMA CLUSTER PRESENCE ACROSS TUMORS
    
    print('\nSupplemental Figure 3c')
    
    #get df from saved csv
    df = pd.read_csv(cwd+'/dfCreated/dfClustCountsN60k7_avg.csv', index_col=0)
    df = df[df.columns[['Perc' in col for col in list(df.columns)]]]
    
    #rename index for swapping P/Rs in chronological order
    xDictRename = {'pre':'A','pos':'B'}
    xAxis = []
    
    #rename xAxis variables to 1P, 1R, etc    
    for x in df.index:
        xLabel = x[2]+xDictRename[x[3:6]]
        xAxis.append(xLabel) #add label to list
    df.index = xAxis
    df = df.sort_index()
    
    #rename index to get P/R
    xDictRename = {'A':'P','B':'R'}
    xAxis2 = []
    
    #rename xAxis variables to 1P, 1R, etc    
    for x in df.index:
        xLabel = x[0]+xDictRename[x[1]]
        xAxis2.append(xLabel) #add label to list
    
    traces = []
    for col in list(df.columns):
        trace = go.Bar(x=xAxis2,y=df[col],name=colDict[col])
        traces.append(trace)
    
    layout = {'xaxis':{'title':'Tumor'},'yaxis':{'title':'Fraction Present'},'barmode':'stack'}
    fig = {'data':traces,'layout':layout}

    pio.write_image(fig,cwd+'/figures/supp3c.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Supplemental Figure 3c saved to 'figures' folder.")


    ##SUPPLEMENTAL FIGURE 3d - ASMA CLUSTER PRESENCE ACROSS PATIENTS
    
    print('\nSupplemental Figure 3d')

    #get df from saved csv
    df = pd.read_csv(cwd+'/dfCreated/dfClustCountsN60k7_avg.csv', index_col=0)
    df = df[df.columns[['Perc' in col for col in list(df.columns)]]]

    df.index = df.index.str.slice(2,3)
    dfPercPt = df.groupby(['patient']).mean()
       
    xAxis = list(dfPercPt.index)
    traces = []
    
    for col in list(dfPercPt.columns):
        trace = go.Bar(x=xAxis,y=dfPercPt[col],name=colDict[col])
        traces.append(trace)
    
    layout = {'xaxis':{'title':'Patient'},'yaxis':{'title':'Fraction Present'},'barmode':'stack'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/supp3d.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Supplemental Figure 3d saved to 'figures' folder.")


    ##SUPPLEMENTAL FIGURE 3e - ASMA CLUSTER PRESENCE ACROSS ANATOMIC SITE
    
    print('\nSupplemental Figure 3e')

    #get df from saved csv
    df = pd.read_csv(cwd+'/dfCreated/dfClustCountsN60k7_all.csv', index_col=0)
    df = df[['pre' in row for row in df.index]]
    
    #generate column list to cluster on based on if there is a % in the column name
    df = df[df.columns[['Perc' in col for col in list(df.columns)]]]
    
    #get clinical csv with stage, anatomic site
    dfClin = pd.read_csv('/Users/blise/Documents/Research/HTAN project/data/grace4/clinicalData.csv',index_col=0)
    dfClin.set_index('sample',inplace=True)
    df['anatomy'] = dfClin.anatomy
    
    #groupby anatomy
    dfPercAnat = df.groupby(['anatomy']).mean()
    dfPercAnat = dfPercAnat.reindex(['oral cavity','oropharynx','larynx']) #set order
    
    #plot
    xAxis = list(dfPercAnat.index)
    traces = []
    
    for col in list(dfPercAnat.columns):
        trace = go.Bar(x=xAxis,y=dfPercAnat[col],name=colDict[col])
        traces.append(trace)
    
    layout = {'xaxis':{'title':'Anatomic Site'},'yaxis':{'title':'Fraction Present'},'barmode':'stack'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/supp3e.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Supplemental Figure 3e saved to 'figures' folder.")
   

    ##SUPPLEMENTAL FIGURE 3f - ASMA CLUSTER PRESENCE ACROSS TNM STAGE
    
    print('\nSupplemental Figure 3f')

    #get df from saved csv
    df = pd.read_csv(cwd+'/dfCreated/dfClustCountsN60k7_all.csv', index_col=0)
    df = df[['pre' in row for row in df.index]]
    
    #generate column list to cluster on based on if there is a % in the column name
    df = df[df.columns[['Perc' in col for col in list(df.columns)]]]
    
    #get clinical csv with stage, anatomic site
    dfClin = pd.read_csv('/Users/blise/Documents/Research/HTAN project/data/grace4/clinicalData.csv',index_col=0)
    dfClin.set_index('sample',inplace=True)
    df['tnm'] = dfClin.tnm.astype(str)
    
    #groupby tnm
    dfPercTnm = df.groupby(['tnm']).mean()
    
    #plot
    xAxis = list(dfPercTnm.index)
    traces = []
    
    for col in list(dfPercTnm.columns):
        trace = go.Bar(x=xAxis,y=dfPercTnm[col],name=colDict[col])
        traces.append(trace)
    
    layout = {'xaxis':{'title':'Anatomic Site'},'yaxis':{'title':'Fraction Present'},'barmode':'stack'}
    fig = {'data':traces,'layout':layout}
    
    pio.write_image(fig,cwd+'/figures/supp3f.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    print("Supplemental Figure 3f saved to 'figures' folder.")
 


if __name__=="__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    table3()
    supp1()
    supp2()
    supp3() 