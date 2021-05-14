#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code author: Katie E. Blise
Date: March, 2021

This .py file contains all code needed to reproduce the figures in the paper:
"Single-Cell Spatial Proteomics Analyses of Head and Neck Squamous Cell Carcinoma Reveal Tumor Heterogeneity and Immune Architectures Associated with Clinical Outcome"


This .py file contains multiple functions:
    getCsvs() - generates a list of all 47 tumor regions analyzed, in correct order for figure generation
    cellClassRename - generates dictionary with gated cell labels converted to biologically relevant cell phenotypes
    countsToCsv() - create csvs with counts of cell types for all tumor regions, to generate figure 1
    kaplan() - creates Kaplan-Meier survival curve
    calculateMixing() - calculates mixing score for all tumor regions, to generate figure 3
    calculateCold() - calculates which tumor regions are labeled 'cold,' to generage figure 3
    pdCounts() - create csvs with counts of PD-1+ and PD-L1+ cells for all tumor regions, to generate figure 3
    kiCounts() - create csvs with counts of Ki-67+ cells for all tumor regions, to generate figure 3
    makeNeighborhoods() - create csv with cellular neighborhood compositions for each seed cell within distance threshold, to generate figure 4
    elbowMethod() - gerate elbow plot to determine optimal number of clusters for k-means clustering
    clusterNeighborhoods() - create csv with results from k-means clustering on cellular neighborhoods, to generate figure 4
    clusterCountPerROI() - create csvs with counts of cellular neighborhood clusters for all tumor regions, to generate figure 4
    fig1() - generates Figures 1A-E
    fig2() - generates Figures 2A-G
    fig3() - generates Figures 3A-H
    fig4() - generates Figures 4B-E
    

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
    #empty list to store csvs in correct order    
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
        
    #print(csvList)
    #print(len(csvList)) #should be 47 files; one per tumor region
    
    return csvList





def cellClassRename():
    
    '''
    This function returns a dicionary of matched cell class IDs (A/B/C/etc) with their corresponding cell class names (CD8+ T cell/CD4+ T-helper cell/B cell/etc).
    
    Input parameters:
        None        
        
    Output:
        returns nameDict = a dictionary of {cell class IDs:cell class names} 
    '''
    
    nameDict = {'A':'CD8+ T Cell','B':'CD4+ T Helper','C':'B Cell','D':'Macrophage','E':'Other Immune','F':'Other Non-Immune','G':'Noise','H':'Neoplastic Tumor','J':'Granulocyte','K':'CD4+ Regulatory T Cell','N':'⍺SMA+ Mesenchymal','X':'Dendritic Cell'}
       
    return nameDict





def countsToCsv(cellList,norm):
    
    '''
    This function calculates the counts (raw or percentage) of cells present in a given list of csvs. It saves the results to csvs - both for all regions and averaged across regions (for a single tumor).
    
    Input parameters:
        cellList = list of cell classes to count; eg ['A','B','C','H']
        norm = True (to get percentage) or False (to get raw count); enter as a boolean NOT a string
        
    Outputs:
        Saves two csvs to the dfCreated folder; one with counts of all regions and one with counts averaged across regions
    '''
    
    #import packages
    import pandas as pd
    import os
    
    #get csvList
    csvList = getCsvs()
    
    #get current directory
    cwd = os.getcwd()

    
    ##create new df and add clinical data to it
    #get clinical df data
    dfClin = pd.read_csv(cwd+'/data/clinicalData.csv', index_col=0)
    
    #empty lists to store clinical data from dfClin
    dtrList = []
    
    #create a new df to store clinical data and counts
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
   
    
    #set up to store counts in value's list of countDict
    countDict = {}
    
    for cell in cellList:
        countDict[cell] = [] #empty list to hold cell counts for each cell type we care about
    
    
    #loop through csvs and calculate counts
    for file in csvList:
        df = pd.read_csv(cwd+'/data/'+file+'.csv',index_col=0)
    
    
        #only include specified cell types in the df
        df = df[[cell in cellList for cell in df['class']]]

    
        #get counts of all cells in the df
        freq = df['class'].value_counts(normalize = norm) #norm specified by input param (raw or %)
    
        for cell in cellList: #only include cells that we want
            
            if cell in freq: #check if the cell type exists in the file
                countDict[cell].append(freq[cell]) #add the frequency of the cell type to the dictionary value's list
            else:
                countDict[cell].append(0) #otherwise add zero as the count
          
        
    #put counts into dfCount
    for c,countList in countDict.items():
        dfCount[c] = countList
    
        
    #csv naming conventions - raw or perc
    if norm == True:
        name = 'Perc_I' #always assume immune cells are being counted
    
        #naming conventions - (I=immune, T=tumor, O=other); only needed for percentage calculation to know what total is out of
        if 'H' in cellList: #includes tumor
            name = name + 'T'
        if ('N' in cellList) or ('F' in cellList): #includes non-immune and non-tumor aka other
            name = name + 'O'    
            
    elif norm == False:
        name = 'Raw'
        
    
    #save dfCount to csv
    dfCount.to_csv(cwd+'/dfCreated/dfCounts'+name+'_all.csv')
    
    #create new dataframe with averaged numbers for every ROI within one patient
    dfCount.index = dfCount.index.str.slice(0,6) #update patient column with just the first six characters in string; so you can use groupby
    
    #groupby same patient (diff ROI) and take the average of the corresponding rows
    dfCountAvg = dfCount.groupby(['patient']).mean()
    
    #save dfAvg to csv
    dfCountAvg.to_csv(cwd+'/dfCreated/dfCounts'+name+'_avg.csv')
    




def kaplan(file,measure):

    '''
    This function creates a Kaplan-Meier survival curve for a given measure.
    It also calculates the p-value using the log-rank test.
    
    Input parameters:
        file = file name saved in dfCreated folder to with data to generate curve
        measure = column name in file to generate curve with (ex: 'A')
        
    Output:
        Returns plt = survival curve with its corresponding p-value
    '''
    
    
    import pandas as pd
    import os
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
    
    #get dictionary to rename cells - used for subfigures
    nameDict = cellClassRename()
    
    #get current directory - used to save subfigures
    cwd = os.getcwd()
  
    #read df
    df = pd.read_csv(cwd+'/dfCreated/'+file+'.csv', index_col=0)
    
    #filter df to only primary tumors
    dfFilt = df[['pre' in row for row in df.index]]
        
    #median of the desired score
    med = dfFilt.median()[measure]
    
    #if median is zero and nothing is lower than zero, then put all zeros into low group and put everyone else in the high group
    if med == 0:
        dfLow = dfFilt[dfFilt[measure] == med] #low group
        dfHigh = dfFilt[dfFilt[measure] > med] #high group
    
    #split into low and high scores - high group has median score
    else:
        dfLow = dfFilt[dfFilt[measure] < med] #low group
        dfHigh = dfFilt[dfFilt[measure] >= med] #high group
    
    
    #fit KM model
    kmf = KaplanMeierFitter()
    
    #create subplot to hold both curves
    plt.figure()
    ax = plt.subplot(111)
    
    #labels for figure legend and title
    if measure == 'mixScore':
        labelL = 'Compartmentalized'
        labelH = 'Mixed'
        title = 'Mixing Score'
    elif 'clust' in measure: #add 1 to each cluster number to match Fig 4B, start cluster count at 1 rather than 0
        labelL = 'Low Cluster '+str(int(measure[5])+1)
        labelH = 'High Cluster '+str(int(measure[5])+1)
        title = 'Cluster '+str(int(measure[5])+1)
    else:
        labelL = 'Low '+nameDict[measure]
        labelH = 'High '+nameDict[measure]
        title = nameDict[measure]
    
    #dfLow curve
    kmf.fit(dfLow['dtr'], event_observed=[1]*len(dfLow), label=labelL)
    kmf.plot(ax=ax,ci_show=False) #can change to True to see confidence intervals
    
    #dfHigh curve
    kmf.fit(dfHigh['dtr'], event_observed=[1]*len(dfHigh), label=labelH)
    kmf.plot(ax=ax,ci_show=False)
    
    plt.ylim(0, 1) #set y axis to go from zero to one
    plt.margins(0) #make (0,0) start in bottom left corner rather than a little bit off
    
    #run logrank test for p-value calculation
    results = logrank_test(dfLow['dtr'],dfHigh['dtr'],[1]*len(dfLow),[1]*len(dfHigh),alpha=0.95)
    
    #title and axis labels
    plt.title(title)
    plt.ylabel('Progression Free Survival (%)')
    plt.xlabel('Days')
    plt.text(x=1150,y=0.7,s='p = '+str(round(results.p_value,3)))
    
    return plt
    




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
    #get clinical df data
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
            
        #print("Mixing score: ",mixScore)
        
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
    
    immCountDict = {}
    
    for file in csvList:
        df = pd.read_csv(cwd+'/data/'+file+'.csv', index_col=0)
    
        counts = df['class'].value_counts()
        
        #reset for each csv
        count = 0
        
        for i in immList:
            if i in counts:
                count = count + counts[i]    
        #add immune count to dictionary for each roi
        immCountDict[file] = count        
    
    #identify csvs with a count of fewer than 2441 immune cells; these are the 'cold' regions
    excludeList = [] #empty list to store cold regions
    for key,value in immCountDict.items():
        if value < 2441: #threshold determined using Keren et al's threshold (adjusting for image size)
            #print(key)
            excludeList.append(key) #add region to the exclude list
    
    #print('Regions to exclude:',excludeList)        
    
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
    It calculates both raw counts as well as percentage counts: for example (PD-1+ CD8 T cells)/(all CD8 T cells)
    Results for % counts are saved to csvs - one for all regions and one averaged across regions (per tumor).
    
    Input parameters:
        None
        
    Outputs:
        Saves two csvs to the dfCreated folder; one with % PD-1+ and PD-L1+ counts of all regions and one with counts averaged across regions    
    '''
    
    import pandas as pd
    import os
    
    #get csvList
    csvList = getCsvs()
    
    #get current directory
    cwd = os.getcwd()

    ##create new df and add clinical data to it
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
    

    #set up to store pd1/pdl1 raw counts in value's list of countDict
    countDict = {}
    dfRaw = pd.DataFrame() #empty df to store raw counts
    dfRaw['patient'] = csvList
    dfRaw.set_index('patient',inplace=True) #set index labels to be the patient column

    
    #list all cell types to count pd-1 and pd-l1 expression
    cellList = ['A','B','C','D','E','J','K','X','F','N','H']   
    
    #create empty lists to hold cell counts for each cell type/PD-status we care about
    for cell in cellList:
        countDict[cell+'_pd1+'] = [] #PD-1+
        countDict[cell+'_pdl1+'] = [] #PD-L1+
    
    #loop through csvs and calculate counts
    for file in csvList:
        #print(file) #to keep track of progress
    
        #get original csv
        df = pd.read_csv(cwd+'/data/'+file+'.csv',index_col=0)
    
    
        for cell in cellList:
            #make smaller df of just the cell type we want
            dfCell = df[(df['class'] == cell)] #note that if cell type does not exist in the df, zeros will automatically be added as the counts
    
            #only count the cells positive for one (or both) of these markers
            #get counts of PD-1
            pd1Pos = len(dfCell[(dfCell['PD1'] == 1)])
    
            #get counts of PD-L1
            pdl1Pos = len(dfCell[(dfCell['PDL1'] == 1)])
        
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
        
    #lists to store % immune cells positive for pd1 and pdl1; will get added to countDict at the end
    allImmPD1List = []
    allImmPDL1List = []
    
    #get original region's csv to get cell counts from
    for pt in ptList:
        #print(pt) #to keep track of progress
        df = pd.read_csv(cwd+'/data/'+pt+'.csv',index_col=0)
    
        #get array of counts for each cell type in the original ROI
        cellCounts = df['class'].value_counts()
            
        #calculate % of each cell
        for c in colList: #loop through each column of the dfRaw file
            
            if c[0] in cellCounts:
                countDict[c].append((dfRaw.loc[pt,c])/(cellCounts[c[0]])) #raw count of cell type/total of that cell type
    
            else:
                countDict[c].append(0) #if cell type is not present in the original ROI, then add zero percent to the countDict for that cell type    
        
        #compute % out of all immune cells
        #get list of cells that are considered immune cells
        immList = ['A','B','C','D','E','J','K','X']

        totImm = 0
        immPD1 = 0
        immPDL1 = 0
        
        for i in immList:
            if i in cellCounts: #only try to add cell count if it actually exists in the original roi
                totImm = totImm + cellCounts[i]
             
            immPD1 = immPD1 + dfRaw.loc[pt,i+'_pd1+']
            immPDL1 = immPDL1 + dfRaw.loc[pt,i+'_pdl1+']
            
        
        allImmPD1List.append(immPD1/totImm)
        allImmPDL1List.append(immPDL1/totImm)
        
    #add all immune cell % counts to countDict
    countDict['allImm_pd1+'] = allImmPD1List
    countDict['allImm_pdl1+'] = allImmPDL1List
        

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
    It calculates both raw counts as well as percentage counts: for example (Ki-67+ CD8 T cells)/(all CD8 T cells)
    Results for % counts are saved to csvs - one for all regions and one averaged across regions (per tumor).
    
    Input parameters:
        None
        
    Outputs:
        Saves two csvs to the dfCreated folder; one with % Ki-67+ counts of all regions and one with counts averaged across regions    
    '''
    
    import pandas as pd
    import os
    
    #get csvList
    csvList = getCsvs()
    
    #get current directory
    cwd = os.getcwd()

    ##create new df and add clinical data to it
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
 
    #set up to store counts in value's list of countDict
    countDict = {}
    dfRaw = pd.DataFrame() #empty df to store raw counts
    dfRaw['patient'] = csvList
    dfRaw.set_index('patient',inplace=True) #set index labels to be the patient column
  
    #list all cell types to count pd-1 and pd-l1 expression
    cellList = ['A','B','C','D','E','J','K','X','F','N','H']   

    #create empty lists to hold cell counts for each cell type/Ki-67-status we are interested in
    for cell in cellList:
        countDict[cell+'_ki67+'] = [] #Ki-67+
    
    #loop through csvs and calculate counts
    for file in csvList:
    
        #get original csv for each tumor region
        df = pd.read_csv(cwd+'/data/'+file+'.csv',index_col=0)
    
        for cell in cellList: #only include cells that we want

            #make smaller df of just the cell type we want
            dfCell = df[(df['class'] == cell)] #note that if cell type does not exist in the df, zeros will automatically be added as the counts
    
            #get counts of Ki-67+ cells
            ki67Pos = len(dfCell[(dfCell['KI67'] == 1)])
        
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
        countDict[c] = [] #empty list to hold cell counts for each cell type/Ki67-status we care about
    
    #lists to store % immune cells positive for Ki67; will get added to countDict at the end
    allImmKi67List = []
    
    #get original region's csv to get cell counts from
    for pt in ptList:
        df = pd.read_csv(cwd+'/data/'+pt+'.csv',index_col=0)
    
        #get array of counts for each cell type in the original region
        cellCounts = df['class'].value_counts()
    
        #calculate % of each cell
        for c in colList: #loop through each column of the dfRaw file
    
            if c[0] in cellCounts:
                countDict[c].append((dfRaw.loc[pt,c])/(cellCounts[c[0]])) #raw count of cell type/total of that cell type
    
            else:
                countDict[c].append(0) #if cell type is not present in the original ROI, then add zero percent to the countDict for that cell type    
    
        #compute % out of all immune cells
        immList = ['A','B','C','D','E','J','K','X']
    
        totImm = 0
        immKi67 = 0
    
        for i in immList:
            if i in cellCounts: #only try to add cell count if it actually exists in the original roi
                totImm = totImm + cellCounts[i]
    
            immKi67 = immKi67 + dfRaw.loc[pt,i+'_ki67+']    
    
        allImmKi67List.append(immKi67/totImm)
    
    #add all immune cell % counts to countDict
    countDict['allImm_ki67+'] = allImmKi67List
    
    #put counts into dfCount
    for c,countList in countDict.items():
        dfCount[c] = countList
    
    #save dfCount to csv; so when you restart kernel you don't have to rerun analysis on all tumors
    dfCount.to_csv(cwd+'/dfCreated/dfKi67CountsPerc_all.csv')
    
    #create new dataframe with averaged numbers for every ROI within one patient
    dfCount.index = dfCount.index.str.slice(0,6) #update patient column with just the first six characters in string; so you can use groupby
    
    #groupby same patient (diff region) and take the average of the corresponding rows
    dfCountAvg = dfCount.groupby(['patient']).mean()
    
    #save dfCountAvg to csv
    dfCountAvg.to_csv(cwd+'/dfCreated/dfKi67CountsPerc_avg.csv')
    




def makeNeighborhoods(seed,distThresh):
    
    '''
    This function counts the neighboring cell types within X distance of the seed cells and saves these counts to a csv.
    Counts are saved both as a raw count and as a percentage.    
    
    Input parameters:
        seed = seed cell type to create neighborhoods for (ex. 'N')
        distThresh = radius (in pixels; 2 px = 1µm) to draw around seed cell to identify neighbors in
    
    Output:
        Saves one csv to the 'dfCreated/' folder with the following naming convention: 'dfNeighborhoodCluster'+seed+str(distThresh)+'.csv'
        This csv is what is used to cluster neighorhoods from in clusterNeighborhoods() function.
        Each row is one seed cell and the columns contain the number of neighbors it has of each cell type.
    '''
    
    import pandas as pd
    from scipy import spatial
    import os
    
    #get csvList
    csvList = getCsvs()
    
    #get current directory
    cwd = os.getcwd()
    
    #empty list to hold dictionaries to create new rows of dfClust
    allNeighList = []
    
    #loop through each ROI in csvList
    for file in csvList:
        #print('Creating neighborhoods for '+file+'...') #to keep track of progress
    
        #load original region's csv
        df = pd.read_csv(cwd+'/data/'+file+'.csv', index_col=0)
    
        #create filtered dataframe without noise
        filt_df = df[(df['class'] != 'G')] #G=noise
    
        #get all possible class valuinges; for later count
        classOptions = sorted(list(filt_df['class'].unique()))
    
        #generate count variable names for final created df's columns
        classCountNames = []
        for c in classOptions: #loop through all possible neighboring cells
            classCountNames.append('count'+c)

        #create np array of just x,y coordinates
        ptsArray = filt_df[['Location_Center_X','Location_Center_Y']].values
    
        #create kdtree
        tree = spatial.KDTree(ptsArray)
        #print('tree created.')
            
        ##get nearest neighbors of seed cells defined by seed param
        #loop through each cell and get its neighbors if it's the right seed
        for i in range(len(ptsArray)):
    
            classType = filt_df['class'].values[i]
    
            #only check neighbors if the seed cell is the desired classType
            if classType == seed:
                neighs = tree.query_ball_point(ptsArray[i], distThresh) #get neighbors within distance
                if i in neighs:
                    neighs.remove(i) #don't include itself as a neighbor
    
                neighClassList = [] #empty list to store neighboring cells' classes
    
                #loop through each neighbor and get its class; add it to neighClassList
                for j in neighs:
                    neighClass = filt_df['class'].values[j]
                    neighClassList.append(neighClass)
                #print(i,neighs,neighClassList)
    
                #get counts of neighboring cell types
                classCounts = []
    
                for c in classOptions: #loop through all possible neighboring cells
                    count = neighClassList.count(c) #count the specified cell type
                    classCounts.append(count) #add the counts of the specified cell type to classCounts
                #print(classCountNames,classCounts)
    
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
                #print('seedDict:',seedDict)
    
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
    #print('Csv created. Check columns to ensure all expected cells are present as columns.')





def elbowMethod(file,steps):
    
    '''
    This function runs the Elbow Method to determine the optimal number of clusters for k-means clustering.
    Note: it is not required to run this function to actually generate any figures, but it was used to determine k=6 for Figure 4 (used steps=15 to determine this).
    It assumes the makeNeighborhoods() function has already been run and dfNeighborhoodClusterN45.csv has been created.
    
    Input parameters:
        file = name of file to run clustering on excluding the .csv (ex. 'dfNeighborhoodClusterN45')
        steps = max number of clusters (k) to test
     
    Output:
        Saves plot to 'figures' folder    
    '''
    
    import pandas as pd
    from matplotlib import pyplot as plt
    from sklearn.cluster import MiniBatchKMeans #use minibatchkmeans when n > 10,000 samples
    import os
    
    #get current directory
    cwd = os.getcwd()
 
    #read in df and generate array with data
    df = pd.read_csv(cwd+'/dfCreated/'+file+'.csv', index_col=0)
            
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
    plt.savefig(cwd+'/figures/ElbowPlot.png',format='png') #saves figure to 'figures' folder as a png file
    #plt.show(graph)
    plt.close()
    print("Elbow plot saved to 'figures' folder.")
 




def clusterNeighborhoods(file,k):
    
    '''
    This function runs k-means clustering on a given neighborhood clustering csv.
    The results are saved to a new csv.
    
    Input parameters:
        file = name of file to run clustering on excluding the .csv (ex. 'dfNeighborhoodClusterN45')
        k = number of clusters; use elbow method to determine optimal number
        
    Outputs:
        Saves csv to the 'dfCreated/' folder as 'dfNighClustered___k_.csv'    
    '''

    import pandas as pd
    from sklearn.cluster import MiniBatchKMeans
    import os
    
    #get current directory
    cwd = os.getcwd()

    #read csv with neighborhood data    
    df = pd.read_csv(cwd+'/dfCreated/'+file+'.csv', index_col=0)
    
    #get lists from original df to add back later after clustering
    roiList = list(df['file']) #get list of regions in order to add to dfNoNa later
    idxList = list(df['index']) #get list of cell indices in order to add to dfNoNa later
        
    #generate column list to cluster on based on if there is a % in the column name
    colList = list(df.columns[['%' in col for col in list(df.columns)]])
    
    #get only features we want to cluster on
    df = df[colList]
    data = df.values
    
    #k-means clustering of cells with k clusters
    kmeans = MiniBatchKMeans(n_clusters=k, init='k-means++', max_iter=300, n_init=10, random_state=0)
    predict = kmeans.fit_predict(data) #fit model to data and predict index (cluster labels)
    #print('prediction complete.')
    
    #add predicted cluster labels to df as a new column
    df['cluster'] = predict
    
    #add original region to df to check which regions are in each cluster
    df['file'] = roiList
    
    #need to add original indices to this df to check which regions are in each cluster; will also need to use regions ID to pair with cell index (since same index could be had by two cells from diff regions)
    df['index'] = idxList #idxList stores the row value of filt_df.iloc[row,column] command
    
    #save df to a csv
    df.to_csv(cwd+'/dfCreated/dfNeighClustered'+file[21:]+'k'+str(k)+'.csv')





def clusterCountPerROI(name):
    
    '''
    This function takes in a clustered csv and creates two csvs with the counts of each seed cell per ROI.
    One csv includes counts across all ROIs. One csv includes counts averaged across ROIs.
    Each csv has both raw count and percentage counts for each cluster.
    
    Input parameters:
        name = name of file containing clustered neighborhood data to count cluster abundance from, excluding the '.csv' (ex. 'dfNeighClusteredN45k6')
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

    #now create new df from dictCounts data and add dtr
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
    df.to_csv(cwd+'/dfCreated/dfClustCounts'+name[16:]+'_all.csv')
    
    #create new dataframe with averaged numbers for every region within one patient
    dfClustCounts.index = dfClustCounts.index.str.slice(0,6)
    
    #groupby same patient (diff ROI) and take the average of the corresponding rows
    dfAvg = dfClustCounts.groupby(['patient']).mean()
    
    #save dfAvg to csv
    dfAvg.to_csv(cwd+'/dfCreated/dfClustCounts'+name[16:]+'_avg.csv')





def fig1():

    '''
    Heterogeneity across patients and tumor regions.
    This function creates figures 1A-1E. 
    All subplots are saved to the 'figures' folder.
    
    Input parameters:
        None
        
    Outputs:
        Saves csvs to 'dfCreated' folder for files containing cell counts used to generate figures.
        Saves plots to 'figures' folder for figures 1A-1E    
    '''
    
    #import packages
    import plotly.graph_objs as go
    import plotly.io as pio
    import pandas as pd
    import os
    from scipy.stats import entropy
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA



    #get dictionary to rename cells - used for subfigures
    nameDict = cellClassRename()
    
    #get current directory - used to save subfigures
    cwd = os.getcwd()
    
    #predefined color scheme; one distinct color per 11 cell types
    cellColors = {'A':'rgb(127,60,141)','B':'rgb(17,165,121)','K':'rgb(57,105,172)','C':'rgb(242,183,1)','D':'rgb(231,63,116)','J':'rgb(128,186,90)','X':'rgb(230,131,16)','E':'rgb(0,134,149)','N':'rgb(207,28,144)','F':'rgb(249,123,114)','H':'rgb(165,170,153)'}

    #dictionary to rename x axis labels
    xDictRename = {'pre':'P','pos':'R'}

   
    
    ##FIGURE 1A - HETEROGENEITY ACROSS PATIENTS AND TUMOR REGIONS - ALL CELL TYPES

    print('\nFigure 1A')
    print('Calculating cell counts...')

    #Calculate counts of ALL cells
    cellList = ['A','B', 'K', 'C', 'D', 'J', 'X', 'E', 'N', 'F', 'H'] #list of cells to count - ALL cells
    countsToCsv(cellList,False) #raw
    countsToCsv(cellList,True) #perc
    
    print("Cell counts calculated. Csvs saved to 'dfCreated' folder.")

    #read cell count file
    df = pd.read_csv(cwd+'/dfCreated/dfCountsPerc_ITO_all.csv',index_col=0)
    
    #empty list to store x axis labels
    xAxis = []
    
    #rename xAxis variables to 1P, 1R, etc    
    for x in df.index:
        xLabel = x[2]+xDictRename[x[3:6]]
        xAxis.append(xLabel) #add label to list
    
    traces = []
    
    for col in list(df.columns)[1:]:
        trace = go.Box(x=xAxis,y=df[col],name=nameDict[col],boxpoints='all',marker={'color':cellColors[col]},line={'color':cellColors[col]}) #note: multiply y= variable by 100 to get out of 100 (as a %)
        traces.append(trace)
    
    layout = {'title':'Intra-Tumoral Heterogeneity (All Cell Types)','xaxis':{'title':'Tumor'},'yaxis':{'title':'Fraction Present'},'boxmode':'group','plot_bgcolor':'rgba(0,0,0,0)'}
    fig = {'data':traces,'layout':layout}    
    pio.write_image(fig,cwd+'/figures/figure1A.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 1A saved to 'figures' folder.")
    
    

    ##FIGURE 1B - HETEROGENEITY ACROSS PATIENTS AND TUMOR REGIONS - IMMUNE CELLS ONLY
    
    print('\nFigure 1B')
    print('Calculating immune cell counts...')

    #Calculate counts of IMMUNE cells only
    cellList2 = ['A','B', 'K', 'C', 'D', 'J', 'X', 'E'] #list of cells to count - only immune cells
    countsToCsv(cellList2,True)

    print("Cell counts calculated. Csvs saved to 'dfCreated' folder.")

    #read cell count file
    df = pd.read_csv(cwd+'/dfCreated/dfCountsPerc_I_all.csv',index_col=0)
        
    #empty list to store x axis labels
    xAxis = []
    
    #rename xAxis variables to 1P, 1R, etc    
    for x in df.index:
        xLabel = x[2]+xDictRename[x[3:6]]
        xAxis.append(xLabel) #add label to list
    
    traces = []
    
    for col in list(df.columns)[1:]:
        trace = go.Box(x=xAxis,y=df[col],name=nameDict[col],boxpoints='all',marker={'color':cellColors[col]},line={'color':cellColors[col]}) #note: multiply y= variable by 100 to get out of 100 (as a %)
        traces.append(trace)
    
    layout = {'title':'Intra-Tumoral Heterogeneity (Immune Cells Only)','xaxis':{'title':'Tumor'},'yaxis':{'title':'Fraction Present'},'boxmode':'group','plot_bgcolor':'rgba(0,0,0,0)'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure1B.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 1B saved to 'figures' folder.")
    
    
    
    ##FIGURE 1C - KL DIVERGENCE
    
    print('\nFigure 1C')
        
    #get df of cell counts - using the perc composition with ALL cells (immune, tumor, other)
    dfAll = pd.read_csv(cwd+'/dfCreated/dfCountsPerc_ITO_all.csv',index_col=0)
    dfAvg = pd.read_csv(cwd+'/dfCreated/dfCountsPerc_ITO_avg.csv',index_col=0)
    
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
    
    #empty list to store KL divergence
    entListA = [] #for intra-tumoral, compares to single tumor's average
    entListB = [] #for intra-patient, compares to single patient's average (P+R)/2
    entListC = [] #for inter-patient, compares to all primary, or all recurrent depending on which one the region is from
    entListD = [] #for inter-patient, compares to cohort average across ALL pts and tumors, regardless of P or R
    
    #create list of rois (patient samples)
    ptList = list(dfAll.index)
    
    #loop through ROIs
    for pt in ptList:
        #get row values for A-H counts for each ROI (per the pt variable)
        countRoi = dfAll.loc[pt,'A':'H'] 
        
        #get character(s) to match other df's index column
        patA = pt[0:6]
        patB = pt[2]
        
        #get row values for A-H counts for each tumor (averaged, use the pat variables)
        countAvgA = dfAvg.loc[patA,'A':'H']
        countAvgB = dfAvgAvg.loc[patB,'A':'H']
        
        if 'pre' in pt:
            countAvgC = dfPrimAvg['A':'H']
            
        elif 'pos' in pt:
            countAvgC = dfRecAvg['A':'H']
    
        countAvgD = dfCohortAvg['A':'H']
        
        #calculate the KL Divergence using the entropy function (log2) for each avg distribution; sample distribution's divergence from the tumor's, patient's, cohort's average (qk) distribution
        entA = entropy(countRoi,qk=countAvgA,base=2) #intra-tumor
        entB = entropy(countRoi,qk=countAvgB,base=2) #intra-patient
        entC = entropy(countRoi,qk=countAvgC,base=2) #inter-patient (p or r)
        entD = entropy(countRoi,qk=countAvgD,base=2) #inter-patient (all)
        
        #store KL divergence in list
        entListA.append(entA)
        entListB.append(entB)
        entListC.append(entC)
        entListD.append(entD)
        
    #create new df with KL divergence info stored in one column
    dfEnt = pd.DataFrame({'Intra-Tumor':entListA,'Intra-Patient':entListB,'Inter-Patient (P or R only)':entListC,'Inter-Patient (all)':entListD},index=ptList)
    
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
    layout = {'title':'Kullback-Leibler Divergence','xaxis':{'title':'Tumor'},'yaxis':{'title':'KL Divergence'},'boxmode':'group'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure1C.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 1C saved to 'figures' folder.")
   
    
    
    ##FIGURE 1D - HIERARCHICAL CLUSTERING OF TUMOR REGION COMPOSITION
    
    print('\nFigure 1D')
    
    #read df
    dfFilt = pd.read_csv(cwd+'/dfCreated/dfCountsPerc_ITO_all.csv',index_col=0)
    
    
    ##color coding - two columns!
    #add a new column for one color - patient number (1-9)
    dfFilt['patient'] = dfFilt.index.str.slice(2,3)
    
    #add a second new column for second color - P(pre)/R(post) (2)
    dfFilt['color2'] = dfFilt.index.str.slice(3,6)
    
    
    #make color vector dictionary; two in total for two color palettes
    palette1 = dict(zip(dfFilt.patient.unique(),sns.color_palette('tab10',n_colors = 9))) ##9 colors for 9 patients
    palette2 = dict(zip(dfFilt.color2.unique(),sns.color_palette("Set3",len(dfFilt.index)))) #Set3 gives light green/yellow combo - adjust in illustrator
    
    #make groupings
    grouping1 = dfFilt.patient #group off of color1 column; renamed to patient so yaxis title says patient rather than color1
    grouping2 = dfFilt.color2 #group off of color2 column
    
    #map palettes onto the corresponding grouping series; names will be used as labels for color columns
    ptColors1 = pd.Series(grouping1,name='Patient ID').map(palette1)
    ptColors2 = pd.Series(grouping2,name='Primary or Recurrent').map(palette2)
    dfColors = pd.concat([ptColors1,ptColors2], axis=1)    
    
    #get only the desired columns aka drop dtr, outcome
    dfFilt = dfFilt[list(dfFilt.columns)[1:-2]]
    
    #rename columns from cell classes to words
    columns = dfFilt.columns #get columns to rename
    #update columns to use real cell classes rather than the IDs
    newCols = []
    for c in columns:
        newCols.append(nameDict[c])
    
    dfFilt.columns = newCols
    
   #plot the clusttermap with the colors; this does the standard_scale (aka normalize [0,1] for each column)
    graph = sns.clustermap(dfFilt,standard_scale=1,cmap='vlag',row_colors=dfColors,yticklabels=True) #yticklabels = true to show all patient labels; row_colors are values of ptColors
    
    #add x-axis label and title to graph
    ax = graph.ax_heatmap
    ax.set_xlabel("\nCell Phenotype")
    ax.set_ylabel("Tumor Region")
    ax.set_title("Hierarchical Clustering of Tumor Region Composition\n\n\n\n\n\n\n\n")
    
    #adjust y limit to fix bug in matplotlib code: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    bottom,top = ax.get_ylim()
    ax.set_ylim(bottom+.5,top-.5)
    
    plt.savefig(cwd+'/figures/figure1D.png',format='png',bbox_inches='tight') #saves figure to 'figures' folder as a png file
    #plt.show(graph)
    plt.close()
    
    print("Figure 1D saved to 'figures' folder.")
    
    
    
    ##FIGURE E - PRINCIPAL COMPONENT ANALYSIS
    
    print('\nFigure 1E')
    
    #read dataframe
    df = pd.read_csv(cwd+'/dfCreated/dfCountsPerc_ITO_all.csv',index_col=0)
    
    #set color palette so it matches other plots; patient colors
    palette = dict(zip(df.index.str.slice(2,3).unique(),sns.color_palette('tab10',n_colors = 9))) ##9 colors for 9 patients
    
    #get roiList, for later dataframe creation, to map to
    roiList = list(df.index)
        
    #separate out features; drop dtr and dtd columns
    feat = df.iloc[:,1:].values
    
    #standardize the features
    feat = StandardScaler().fit_transform(feat)
    
    
    ##DETERMINE THE TARGET; could be pre vs post or patient id, etc; comment/uncomment the appropriate lines of code for the desired target
    
    #target = patient ID
    df['target'] = [roi[2] for roi in df.index.values]
        
    
    #do PCA
    pca = PCA(n_components = 2) #generate model based on top two PCs
    principalComponents = pca.fit_transform(feat) #fit and transform features
    
    #create new dataframe holding the principal components; need to specify index so that you can map the target column in next step
    dfPrinc = pd.DataFrame(data=principalComponents, columns=['PC1','PC2'], index = roiList)
    
    #concatenate the principal component df and the target column (this is what gets color coded in final plot)
    dfFinal = pd.concat([dfPrinc,df[['target']]], axis=1)
    
    #rename target column name to appropriate name for plotting - to patient
    dfFinal = dfFinal.rename(columns={'target':'Patient'})
    
    
    #VISUALIZE using seaborn, color according to the target column
    sns.set(rc={'figure.figsize':(6,6)})
    graph = sns.scatterplot(dfFinal['PC1'], dfFinal['PC2'], hue=dfFinal['Patient'], legend='full',palette=palette)
    graph.legend(loc='lower right',bbox_to_anchor=(1.25, 0.5),ncol=1)
    graph.set_title("Principle Component Analysis")
    
    plt.savefig(cwd+'/figures/figure1E.png',format='png',bbox_inches='tight') #saves figure to 'figures' folder as a png file
    #plt.show(graph)
    plt.close()    
    
    print("Figure 1E saved to 'figures' folder.")
        
        
        
    
            
def fig2():
    
    '''
    Tumor cellular composition changes following therapy and composition associated with survival.
    This function creates figures 2A-2G. 
    All subplots are saved to the 'figures' folder.
    
    Note: This function assumes countsToCsvs() has already been run (this happens in fig1() function), as it relies on count csvs to generate figures.
    
    Input parameters:
        None
        
    Outputs:
        Saves plots to 'figures' folder for figures 2A-2G
    '''

    #import packages
    import pandas as pd
    import os
    import plotly.graph_objs as go
    import plotly.io as pio
    from scipy.stats import ttest_rel, wilcoxon
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA



    #get dictionary to rename cells - used for subfigures
    nameDict = cellClassRename()
    
    #get current directory - used to save subfigures
    cwd = os.getcwd()



    ##FIGURE 2A - PRIMARY VS RECURRENT CELL COUNTS
    
    print('\nFigure 2A')
    
    #read in counts df
    df = pd.read_csv(cwd+'/dfCreated/dfCountsPerc_ITO_avg.csv',index_col=0)
    
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
    
    
    layout = {'title':'Primary vs Recurrent TiME Cell Counts' ,'xaxis':{'title':'Cell Phenotype'},'yaxis':{'title':'% Composition (Averaged Across Regions)'},'boxmode':'group'}    
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure2A.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Figure 2A saved to 'figures' folder.")
    
    #check for significance between pre/post tumors and print the p-value, paired tests
    print('Calculating p-values...')
    for col in dfPre.columns: #for each column in the df
    
        valueList = [] #overall list to store two lists to compare (pre and post); reset with each new column
    
        for df in dfList: #for each df (pre and then post)
            valList = list(df[col])
            valueList.append(valList) #add pre or post's value list to overall valueList
        
        print('\n',nameDict[col]+':')
        
        #non-parametric paired t-test (for B cell counts, whose difference does NOT follow a normal distribution)
        if col == 'C':
            print('Non-parametric (Wilcoxon signed rank) paired t-test (two-tailed) p-value:',wilcoxon(valueList[0],valueList[1])[1])
        
        #all other cell types' differences DO follow a normal distribution, so use parametric paired t-test
        else:
            print('Parametric paired t-test (two-tailed) p-value:',ttest_rel(valueList[0],valueList[1])[1])



    ##FIGURE 2B - PRIMARY VS RECURRENT CELL COUNTS

    print('\nFigure 2B')

    #read in counts df
    df = pd.read_csv(cwd+'/dfCreated/dfCountsPerc_ITO_avg.csv',index_col=0)
    
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
    fig, ax =plt.subplots(2,6,figsize=(25,25))
    
    #range through 12 to set variables for placement of graph, where in the df to get data from
    for n in range(12):
        
        #set placement of graph
        #row
        if n < 6:
            a = 0
        else:
            a = 1
        #column
        if n < 6:
            b = n
        else:
            b = n-6
        
        #set iloc to draw data from
        i = n*2
        j = i+1
        
        if i < 22:
            
            ax[a,b].set_title(titleList[n],fontdict={'fontsize':14}) #set title of each subplot
            ax[a,b].set_ylabel("% Composition (Averaged Across Regions)") #set y axis label of each subplot        
            #make lineplot; use tab10 color scheme to match before
            sns.lineplot(data=dfFinal.iloc[[i,j],:]*100,dashes=False,palette='tab10',ax=ax[a,b],legend=False,marker="o")
    
        #skip the final plot that does not exist 
        else:
            break
    
        
    #final adjustments
    fig.delaxes(ax= ax[1,5]) #delete the 12th ax plot since there are only 11 cells
    fig.suptitle("Primary vs Recurrent TiME Cell Counts", fontsize=20,va='bottom') #add global title to figure
    fig.subplots_adjust(top=.95) #make title lower
    ax[1,4].legend(['1', '2', '3','4','5','6','7','8','9'],bbox_to_anchor=(1.9,0.75),ncol=1,title='Patient',prop={'size': 18}) #create legend for the final ax 
    
    plt.savefig(cwd+'/figures/figure2B.png',format='png') #saves figure to 'figures' folder as a png file
    plt.close()
    print("Figure 2B saved to 'figures' folder.")



    ##FIGURE 2C - PRINCIPAL COMPONENT ANALYSIS
    
    print('\nFigure 2C')
    
    #read dataframe
    df = pd.read_csv(cwd+'/dfCreated/dfCountsPerc_ITO_all.csv',index_col=0)
    
    #get roiList, for later dataframe creation, to map to
    roiList = list(df.index)
    
    #separate out features; drop dtr and dtd columns
    feat = df.iloc[:,1:].values
    
    #standardize the features
    feat = StandardScaler().fit_transform(feat)
    
    #target = pre vs post for _all regions csv
    df['target'] = [roi[3:6] for roi in df.index.values]
    
    #do PCA
    pca = PCA(n_components = 2) #make model based on top two PCs
    principalComponents = pca.fit_transform(feat) #fit and transform features
    
    #create new dataframe holding the principal components; need to specify index so that you can map the target column in next step
    dfPrinc = pd.DataFrame(data=principalComponents, columns=['PC1','PC2'], index = roiList)
    
    #concatenate the principal component df and the target column (this is what gets color coded in final plot)
    dfFinal = pd.concat([dfPrinc,df[['target']]], axis=1)
    
    #rename target column name to appropriate name for plotting - to patient
    dfFinal = dfFinal.rename(columns={'target':'Tumor'})
    
    #empty list to store x axis labels
    ptList = []
    
    #rename xAxis variables to 1P, 1R, etc
    xDictRename = {'pre':'Primary','pos':'Recurrent'}
    
    for x in dfFinal.index:
        xLabel = xDictRename[x[3:6]]
        ptList.append(xLabel) #add label to list
    dfFinal.index = ptList
    dfFinal.Tumor = ptList
    
    #VISUALIZE using seaborn, color according to the target column
    sns.set(rc={'figure.figsize':(6,6)})
    graph = sns.scatterplot(dfFinal['PC1'], dfFinal['PC2'], hue=dfFinal['Tumor'], legend='full')#,palette=palette)
    graph.set_title("Principle Component Analysis")
    
    plt.savefig(cwd+'/figures/figure2C.png',format='png') #saves figure to 'figures' folder as a png file
    #plt.show(graph)    
    plt.close()
    print("Figure 2C saved to 'figures' folder.")        
    
    
    
    ##FIGURE 2D - HIERARCHICAL CLUSTERING OF AVERAGE CHANGE IN TIME CELLULAR COMPOSITION
    
    print('\nFigure 2D')
    
    #read df
    df = pd.read_csv(cwd+'/dfCreated/dfCountsPerc_ITO_avg.csv',index_col=0)
    
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
    dfClin2 = dfClin.copy()
    dfClin2['sample'] = dfClin2['sample'].str.slice(2,3)
    dfClin2 = dfClin2.drop_duplicates()
    dfClin2.set_index('sample',inplace=True)
    
    dfDiff['patient'] = dfDiff.index
    dfDiff['tx'] = dfClin2['tx']
    
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
    
    #plot the clusttermap with the colors; this does standard_scale ([0,1] normalization)
    graph = sns.clustermap(dfDiff,standard_scale=1,cmap='vlag',row_colors=dfColors,yticklabels=True)
    
    #add x-axis label and title to graph
    ax = graph.ax_heatmap
    ax.set_xlabel("\nCell Phenotype")
    ax.set_ylabel("Patient")
    ax.set_title("Hierarchical Clustering of Average Change in TiME Cellular Composition\n\n\n\n\n\n\n\n")
    
    #adjust y limit to fix bug in matplotlib code
    bottom,top = ax.get_ylim()
    ax.set_ylim(bottom+.5,top-.5)
    
    plt.savefig(cwd+'/figures/figure2D.png',format='png',bbox_inches='tight') #saves figure to 'figures' folder as a png file
    #plt.show(graph)
    plt.close()
    print("Figure 2D saved to 'figures' folder.")
                
    
    
    ##FIGURE 2E - SURVIVAL CURVE, CD8+ T CELL
    
    print('\nFigure 2E')
    
    #run kaplan function to generate survival curve and return plot to save
    plt = kaplan('dfCountsPerc_ITO_avg','A') #A = CD8+ T cell
    
    plt.savefig(cwd+'/figures/figure2E.png',format='png') #saves figure to 'figures' folder as a png file
    #plt.show()
    plt.close()
    print("Figure 2E saved to 'figures' folder.")        

    
    
    ##FIGURE 2F - SURVIVAL CURVE, MACROPHAGE
    
    print('\nFigure 2F')
    
    #run kaplan function to generate survival curve and return plot to save
    plt = kaplan('dfCountsPerc_ITO_avg','D') #D = macrophage
    
    plt.savefig(cwd+'/figures/figure2F.png',format='png') #saves figure to 'figures' folder as a png file
    #plt.show()
    plt.close()
    print("Figure 2F saved to 'figures' folder.")        

    

    ##FIGURE 2G - SURVIVAL CURVE, OTHER NON-IMMUNE
    
    print('\nFigure 2G')
    
    #run kaplan function to generate survival curve and return plot to save
    plt = kaplan('dfCountsPerc_ITO_avg','F') #F = other non-immune
    
    plt.savefig(cwd+'/figures/figure2G.png',format='png') #saves figure to 'figures' folder as a png file
    #plt.show()
    plt.close()
    print("Figure 2G saved to 'figures' folder.")        




    
def fig3():
    
    '''
    Mixing score quantifies the spatial organization of tumors.
    This function creates figures 3A-3H. 
    All subplots are saved to the 'figures' folder.
    
    Note: This function assumes countsToCsvs() has already been run (this happens in fig1() function), as it relies on count csvs to generate figures 3C and 3D.
    
    Input parameters:
        None
        
    Outputs:
        Saves csvs to 'dfCreated' folder for files containing mixing scores, PD-1+ counts, PD-L1+ counts, and Ki-67+ counts used to generate figures.
        Saves plots to 'figures' folder for figures 3A-3H    
    '''
  
    #import packages    
    import pandas as pd
    import os
    import plotly.graph_objs as go
    import plotly.io as pio
    from scipy.stats import mannwhitneyu, ttest_ind
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    import seaborn as sns
    import matplotlib.pyplot as plt



    #get dictionary to rename cells - used for subfigures
    nameDict = cellClassRename()
    
    #get current directory - used to save subfigures
    cwd = os.getcwd()
    
    #predefined color scheme; one distinct color per 11 cell types
    cellColors = {'A':'rgb(127,60,141)','B':'rgb(17,165,121)','K':'rgb(57,105,172)','C':'rgb(242,183,1)','D':'rgb(231,63,116)','J':'rgb(128,186,90)','X':'rgb(230,131,16)','E':'rgb(0,134,149)','N':'rgb(207,28,144)','F':'rgb(249,123,114)','H':'rgb(165,170,153)'}
    
    #set color dict to map colors from for plotting spatial organization
    colorDict = {'Cold':'#E377C2','Mixed':'#FF7F03','Compartmentalized':'#2CA02C'}

    #rename xAxis variables to 1P, 1R, etc
    xDictRename = {'pre':'P','pos':'R'}

    #create function to map spatial labels into dfs within figure 3 subplots
    def map_mix(mix):
    
        if mix < 0.1034: #value corresponds to median of primary tumor regions
            return 'Compartmentalized'
        elif mix > 0.1034:
            return 'Mixed'



    ##FIGURE 3A - VISUALIZE MIXED AND COMPARTMENTALIZED TUMOR REGIONS
    
    print('\nFigure 3A')
    
    #example regions: high mixing - pt7post_roi2, low mixing = pt8pre_roi3
    files = ['pt7post_roi2','pt8pre_roi3']
    
    titleDict = {'pt7post_roi2':'Mixed','pt8pre_roi3':'Compartmentalized'}
    
    for file in files:
            
        #read original csv file
        df = pd.read_csv(cwd+'/data/'+file+'.csv', index_col=0)
    
        #list of lists - cell types to plot - immune only (I), tumor only (T), immune&tumor (IT)
        classListBig = [['E','A','B','K','C','D','J','X'],['H'],['H','E','A','B','K','C','D','J','X']]
    
        classListBigDict = {0:'I',1:'T',2:'IT'} #dict for figure file naming based on i in next loop
    
        #loop through each cell class list to plot
        for i in range(len(classListBig)):    
    
            classList = classListBig[i]
    
            #create filtered dataframe with specific cell types present
            dfFilt = df[df['class'].isin(classList)]
    
            #plot cells in each of the classLists
            traces = []
            for celltype in classList: #for each cell type
                filtered_df = dfFilt[dfFilt['class'] == celltype] #filtered_df is a dataframe of each specific cell type
                #create one trace per cell type using Location_Center_X/Y as (x,y) coordinates
                trace = go.Scatter(x = filtered_df.Location_Center_X,y = filtered_df.Location_Center_Y,mode = 'markers',name = nameDict[celltype],marker={'color':cellColors[celltype]})
                traces.append(trace)
    
            #generate plot
            layout = {'title':titleDict[file],'showlegend':True,'autosize':False,'width':800,'height':800,'xaxis':{'visible':False},'yaxis':{'visible':False}}
            fig = {'data':traces,'layout':layout}
            pio.write_image(fig,cwd+'/figures/figure3A_'+titleDict[file]+'_'+classListBigDict[i]+'.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 3A saved to 'figures' folder.")

    
    
    ##FIGURE 3B - MIXING SCORE ACROSS TUMOR REGIONS
    
    print('\nFigure 3B')
    
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
    coldList = ['pt2post_roi1','pt4pre_roi1','pt4pre_roi2','pt9post_roi3'] #this is a pre-defined list based on # of immune cells (fewer than 2441 imm)
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
    
    #trace with all data points, box
    traceData = go.Box(x=xAxis,y=df['mixScore'],boxmean=True,boxpoints='all',name='data',jitter=1)
    
    #create traces to store all trace variables
    traces = [traceData]
    
    #create new trace of markers only for color-coding by spatial org
    for cat in df['spatial'].unique():
    
        #get a sub-df with rois that all have the same spatial category
        dfSpat = df[df['spatial'] == cat]
    
        #create marker trace for that spatial category and add it to traces list
        traceSpat = go.Scatter(x=list(dfSpat.index),y=dfSpat['mixScore'],mode='markers',marker = {'size':10,'color':colorDict[cat]},name=cat)
        traces.append(traceSpat)
    
    #create annotations with overall averaged mixing score label
    #get overall mixing score designation from average csv 
    file2 = 'dfMixingScoresExcludeCold_avg'
    dfAvg = pd.read_csv(cwd+'/dfCreated/'+file2+'.csv', index_col=0)
    #add spatial column
    dfAvg['spatial'] = dfAvg['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold tumors where the regions are all cold
    coldList = ['pt2post','pt4pre'] #this is a pre-defined list based on # of immune cells (fewer than 2441 imm)
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
    layout = {'title':'Mixing Score Across Tumor Regions','xaxis':{'title':'Tumor'},'yaxis':{'title':'Mixing Score'},'annotations':annotations}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure3B.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Figure 3B saved to 'figures' folder.")
    
    
    
    ##FIGURE 3C - NEOPLASTIC TUMOR CELL COUNT COLORED BY MIXING SCORE
    
    print('\nFigure 3C')
    
    #read csv with counts in it - this csv gets created in fig1() function
    dfCounts = pd.read_csv(cwd+'/dfCreated/dfCountsRaw_all.csv', index_col=0)
    
    #H = neoplastic tumor cells
    count = dfCounts['H']
    
    #sort regions by count in ascending order
    count = count.sort_values()
    
    #get labels of the regions; x axis
    roiLabels = list(count.index)
    
    #get counts to plot; y axis
    countValues = list(count.values)
    
    #read csv with mixing scores
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv',index_col=0)
    
    #get just the allImmune mix score column
    mix = dfMix['mixScore']
    
    #empty list to store colors according to mixing score
    colorList = []
    mixLabel = []
    
    #loop through each region and assign color based on its mixing score and immune count (for cold only)
    for roi in roiLabels:
    
        #if not cold; 2441 cells determined by using same % of immune cells in µm^2 as Keren et al. paper
        if roi != 'pt2post_roi1' and roi != 'pt4pre_roi1' and roi != 'pt4pre_roi2' and roi != 'pt9post_roi3':
    
            #compartmentalized (median of averaged mixing scores)
            if mix[roi] < 0.1034:
                colorList.append('#2CA02C')
                mixLabel.append('Compartmentalized')
            #mixed
            else:
                colorList.append('#FF7F03')
                mixLabel.append('Mixed')
    
        #cold; fewer than 2441 cells
        else:
            colorList.append('#E377C2')
            mixLabel.append('Cold')
    
    #make figure
    trace = [go.Bar(x=roiLabels,y=countValues,marker={'color':colorList})]
    
    #make annotations for mixing labels to appear on plot
    annotC = dict(x=38,y=39000,text='Compartmentalized',showarrow=False,font={'color':'#2CA02C'})
    annotM = dict(x=38,y=37000,text='Mixed',showarrow=False,font={'color':'#FF7F03'})
    annotC2 = dict(x=38,y=35000,text='Cold',showarrow=False,font={'color':'#E377C2'})
    annotList = [annotC,annotC2,annotM]
    
    layout = {'title':'Neoplastic Tumor Cell Count Colored by Mixing Score Per Region','xaxis':{'title':'Region','automargin':True},'yaxis':{'title':'Neoplastic Tumor Cell Count'},'annotations':annotList}
    fig = {'data':trace,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure3C.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 3C saved to 'figures' folder.")
    
    
    
    ##FIGURE 3D - ASMA MESENCHYMAL CELL COUNT COLORED BY MIXING SCORE
    
    print('\nFigure 3D')
    
    #read csv with counts in it - this csv gets created in fig1() function
    dfCounts = pd.read_csv(cwd+'/dfCreated/dfCountsRaw_all.csv', index_col=0)
    
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
    
    #get just the allImmune mix score column
    mix = dfMix['mixScore']
    
    #empty list to store colors according to mixing score
    colorList = []
    mixLabel = []
    
    #loop through each region and assign color based on its mixing score and immune count (for cold only)
    for roi in roiLabels:
    
        #if not cold; 2441 cells determined by using same % of immune cells in µm^2 as Keren et al. paper
        if roi != 'pt2post_roi1' and roi != 'pt4pre_roi1' and roi != 'pt4pre_roi2' and roi != 'pt9post_roi3': #determined by calculateCold() function
    
            #compartmentalized (median of averaged mixing scores)
            if mix[roi] < 0.1034:
                colorList.append('#2CA02C')
                mixLabel.append('Compartmentalized')
            #mixed
            else:
                colorList.append('#FF7F03')
                mixLabel.append('Mixed')
    
        #cold; fewer than 2441 cells
        else:
            colorList.append('#E377C2')
            mixLabel.append('Cold')
    
    #make figure
    trace = [go.Bar(x=roiLabels,y=countValues,marker={'color':colorList})]
    
    #make annotations for mixing labels to appear on plot
    annotC = dict(x=37,y=3000,text='Compartmentalized',showarrow=False,font={'color':'#2CA02C'})
    annotM = dict(x=37,y=2850,text='Mixed',showarrow=False,font={'color':'#FF7F03'})
    annotC2 = dict(x=37,y=2700,text='Cold',showarrow=False,font={'color':'#E377C2'})
    annotList = [annotC,annotC2,annotM]
    
    layout = {'title':'⍺SMA+ Mesenchymal Cell Count Colored by Mixing Score Per Region','xaxis':{'title':'Region','automargin':True},'yaxis':{'title':'⍺SMA+ Mesenchymal Cell Count'},'annotations':annotList}
    fig = {'data':trace,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure3D.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 3D saved to 'figures' folder.")
    
    
    
    ##FIGURE 3E - PRINCIPAL COMPONENT ANALYSIS
    
    print('\nFigure 3E')
    
    #set color palette so it matches other plots; mixing score colors
    palette = dict(Compartmentalized = '#2CA02C',Mixed = '#FF7F03',Cold = '#E377C2')
    
    #read counts csv
    df = pd.read_csv(cwd+'/dfCreated/dfCountsPerc_ITO_all.csv',index_col=0)
    
    #get roiList, for later dataframe creation, to map to
    roiList = list(df.index)
    
    #separate out features; drop dtr column
    feat = df.iloc[:,1:].values
    
    #standardize the features
    feat = StandardScaler().fit_transform(feat)
        
    #add a column on mixing df for mixed, compartmentalized, cold - to see intra-tumor heterogeneity    
    #read in mixing score csv - all regions
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix)) #map_mix defined above
    
    #adjust df to also account for cold ROIs - all ROIs
    coldList = ['pt2post_roi1','pt4pre_roi1','pt4pre_roi2','pt9post_roi3'] #this is a pre-defined list based on # of immune cells (fewer than 2441 imm)
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
    
    #VISUALIZE using seaborn, color according to the target column
    sns.set(rc={'figure.figsize':(7,7)})
    graph = sns.scatterplot(dfFinal['PC1'], dfFinal['PC2'], hue=dfFinal['Spatial Organization'], legend='full',palette=palette)
    graph.set_title("Principle Component Analysis")
    
    plt.savefig(cwd+'/figures/figure3E.png',format='png') #saves figure to 'figures' folder as a png file
    #plt.show(graph)
    plt.close()
    print("Figure 3E saved to 'figures' folder.")    
    
    
    ##FIGURE 3F - NEOPLASTIC TUMOR PD-1 AND PD-L1 EXPRESSION BY SPATIAL ARCHITECTURE
    
    print('\nFigure 3F') 
    
    #calculate pd-1, pd-l1 positive cell counts
    print('Calculating PD-1+ and PD-L1+ cell counts for each tumor region...')
    pdCounts()
    print("PD-1+ and PD-L1+ cell counts calculated. Csvs saved to 'dfCreated' folder.")
    
    #read pd1/pdl1 counts csv
    df = pd.read_csv(cwd+'/dfCreated/dfPD1PDL1CountsPerc_all.csv', index_col=0)
    
    #only want tumor counts for fig 3F - aka H_pd1+ and H_pdl1+; can adjust this to get other cell counts
    dfCount = df.loc[:,'H_pd1+':'H_pdl1+'].copy() #copy to avoid chaining
    
    #read in mixing score df for all regions
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix)) #map_mix defined above
    
    #adjust df to also account for cold regions
    coldList = ['pt2post_roi1','pt4pre_roi1','pt4pre_roi2','pt9post_roi3'] #based on # of immune cells (fewer than 2441 imm)
    for roi in coldList:
            dfMix.loc[roi,'spatial'] = 'Cold'
    
    #now add the spatial header into the target column in the dfCount
    dfCount['spatial'] = list(dfMix['spatial'])
    
    #separate out mixed and compartmentalized based on spatial column
    dfMixed = dfCount[['Mixed' in row for row in dfCount['spatial']]]
    dfCompart = dfCount[['Compartmentalized' in row for row in dfCount['spatial']]]
    
    #get values to plot for compartmentalized trace
    xCom = []
    for i in range(len(dfCompart)):
        xCom.append('PD-1')
    for i in range(len(dfCompart)):
        xCom.append('PD-L1')
    
    yComPD1 = list(dfCompart.loc[:,'H_pd1+'].values*100) #values are out of 100%
    yComPDL1 = list(dfCompart.loc[:,'H_pdl1+'].values*100)
    yCom = yComPD1 + yComPDL1
    
    #get values to plot for mixed trace
    xMix = []
    for i in range(len(dfMixed)):
        xMix.append('PD-1')
    for i in range(len(dfMixed)):
        xMix.append('PD-L1')
    
    yMixPD1 = list(dfMixed.loc[:,'H_pd1+'].values*100)
    yMixPDL1 = list(dfMixed.loc[:,'H_pdl1+'].values*100)
    yMix = yMixPD1 + yMixPDL1   
    
    traceC = go.Box(x=xCom,y=yCom,name='Compartmentalized',boxpoints='all',marker={'color':'#2CA02C'},line={'color':'#2CA02C'})
    traceM = go.Box(x=xMix,y=yMix,name='Mixed',boxpoints='all',marker={'color':'#FF7F03'},line={'color':'#FF7F03'})
    traces = [traceC,traceM]
    
    layout = {'title':'Neoplastic Tumor PD-1 and PD-L1 Expression By Spatial Architecture','xaxis':{'title':'Immunoregulatory Protein'},'yaxis':{'title':'% of Neoplastic Cells Positive'},'boxmode':'group'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure3F.png',format='png',scale=2) #saves figure to 'figures' folder as a png file
    
    print("Figure 3F saved to 'figures' folder.")
   
    #calculate p-values
    print('Calculating statistical significance...')
    print('PD-1 Mann-Whitney U test (one-sided, p-value):',mannwhitneyu(dfMixed.loc[:,'H_pd1+'],dfCompart.loc[:,'H_pd1+'],alternative='greater')[1]) #non-parametric independent t-test
    print('PD-L1 Mann-Whitney U test (one-sided, p-value):',mannwhitneyu(dfMixed.loc[:,'H_pdl1+'],dfCompart.loc[:,'H_pdl1+'],alternative='greater')[1]) #non-parametric independent t-test
              
        
    
    ##FIGURE 3G - LYMPHOCYTE KI-67 EXPRESSION BY SPATIAL ARCHITECTURE
    
    print('\nFigure 3G') 
    
    #calculate KI-67 positive cell counts
    print('Calculating Ki-67+ cell counts for each tumor region...')
    kiCounts()
    print("Ki-67+ cell counts calculated. Csvs saved to 'dfCreated' folder.")

    #read ki-76 counts csv
    df = pd.read_csv(cwd+'/dfCreated/dfKi67CountsPerc_all.csv', index_col=0)
    
    #only keep A, B, C, and K columns = lymphocytes only
    cellsToKeep = ['A','B','C','K']
    dfCount = df[[col for col in list(df.columns) if col[0] in cellsToKeep]]
    
    #sum across column series
    dfKiSum = pd.DataFrame(dfCount.sum(axis=1),columns=['Ki-67+'])
     
    #read in appropriate mixing score df (either all or avg)
    dfMix = pd.read_csv(cwd+'/dfCreated/dfMixingScores_all.csv', index_col=0)
    
    dfMix['spatial'] = dfMix['mixScore'].apply(lambda mix:map_mix(mix))
    
    #adjust df to also account for cold ROIs - all ROIs
    coldList = ['pt2post_roi1','pt4pre_roi1','pt4pre_roi2','pt9post_roi3'] #this is a pre-defined list based on # of immune cells (fewer than 2441 imm)
    for roi in coldList:
            dfMix.loc[roi,'spatial'] = 'Cold'
    
    #now add the spatial header into the spatial column to the Kisum df
    dfKiSum['spatial'] = list(dfMix['spatial'])
    
    #separate out mixed and compartmentalized based on spatial column
    dfMixed = dfKiSum[['Mixed' in row for row in dfKiSum['spatial']]]
    dfCompart = dfKiSum[['Compartmentalized' in row for row in dfKiSum['spatial']]]
    
    traceM = go.Box(y=dfMixed.iloc[:,0]*100,name='Mixed',boxpoints='all',marker={'color':'#FF7F03'},line={'color':'#FF7F03'})
    traceC = go.Box(y=dfCompart.iloc[:,0]*100,name='Compartmentalized',boxpoints='all',marker={'color':'#2CA02C'},line={'color':'#2CA02C'})
    
    traces = [traceC,traceM]
    layout = {'title':'Lymphocyte Ki-67 Expression by Spatial Architecture','xaxis':{'title':'Spatial Architecture'},'yaxis':{'title':'% of Lymphocytes Positive'},'boxmode':'group'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure3G.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 3G saved to 'figures' folder.")
   
    #calculate p-values
    print('Calculating statistical significance...')
    #Run a two-tailed independent, parametric t-test (since lymphocyte diff follows a normal distribution) with order (compart,mixed)
    #H0: C <= M
    #HA: C > M
    #then divide the p value by 2 to get a one-sided t-test
    #if this p value is less than alpha (0.05), then you can reject the null in favor of the alternative
    pTwo = ttest_ind(dfCompart.iloc[:,0],dfMixed.iloc[:,0],equal_var=False)[1] #returns results from TWO-TAIL independent (parametric) t-test
    pVal = pTwo/2 #must divide p by two to get results of one-tailed test
    print('One-sided t-test, p-value:',pVal)
 
    
    
    ##FIGURE 3H - SURVIVAL CURVE, MIXING SCORE
    
    print('\nFigure 3H') 
    
    #run kaplan function to generate survival curve and return plot to save
    plt = kaplan('dfMixingScoresExcludeCold_avg','mixScore') #F = other non-immune
    
    plt.savefig(cwd+'/figures/figure3H.png',format='png') #saves figure to 'figures' folder as a png file
    #plt.show()
    plt.close()
    print("Figure 3H saved to 'figures' folder.")        

    
    
    
    
def fig4():
    
    '''
    αSMA+ mesenchymal cellular neighborhood clustering.
    This function creates figures 4B-4E. 
    All subplots are saved to the 'figures' folder.
    Note: Figure 4A was created using BioRender.com.
        
    Input parameters:
        None
        
    Outputs:
        Saves csvs to 'dfCreated' folder for files containing cellular neighborhood data used to generate figures.
        Saves plots to 'figures' folder for figures 4B-4E    
    '''
  
    #import packages    
    import pandas as pd
    import os
    import plotly.io as pio
    


    #get dictionary to rename cells - used for subfigures
    nameDict = cellClassRename()
    
    #get current directory - used to save subfigures
    cwd = os.getcwd()
    
    #predefined color scheme; one distinct color per 11 cell types
    cellColors = {'A':'rgb(127,60,141)','B':'rgb(17,165,121)','K':'rgb(57,105,172)','C':'rgb(242,183,1)','D':'rgb(231,63,116)','J':'rgb(128,186,90)','X':'rgb(230,131,16)','E':'rgb(0,134,149)','N':'rgb(207,28,144)','F':'rgb(249,123,114)','H':'rgb(165,170,153)'}



    ##FIGURE 4B - AVERAGE αSMA+ MESENCHYMAL CELL NEIGHBORHOOD CLUSTER COMPOSITION

    print('\nFigure 4B')

    #make neighborhoods for αSMA+ seed cells (N), distThresh = 22.5µm = 45px
    print('Creating αSMA+ cellular neighborhoods...')
    makeNeighborhoods('N',45)
    print("Neighborhoods created. Csv saved to 'dfCreated' folder.")

    #optionally run elbowMethod() function to determine optimal number of clusters
    print('Running Elbow Method to determine optimal number of clusters...')
    elbowMethod('dfNeighborhoodClusterN45',15) #move forward with k=6 based on these results
    
    #cluster αSMA+ cellular neighborhoods based on composition, k=6
    print('Clustering αSMA+ cellular neighborhoods...')
    clusterNeighborhoods('dfNeighborhoodClusterN45',6)
    print("Neighborhoods clustered. Csv saved to 'dfCreated' folder.")
        
    #read csv
    df = pd.read_csv(cwd+'/dfCreated/dfNeighClusteredN45k6.csv', index_col=0)
    
    #filter df to not include the index or file columns
    df = df.drop(['index','file'],axis=1)
    
    #groupby cluster column and take the averages of all of the other columns for each group
    dfCluster = df.groupby(['cluster']).mean()
    
    #transpose df so columns are clusters and rows are cell types
    dfClusterT = dfCluster.T
    
    #rename clusters to be 1-6 rather than 0-5
    dfClusterT = dfClusterT.rename(columns={0:'1',1:'2',2:'3',3:'4',4:'5',5:'6'})
    
    #re-order the rows of dfClusterT such that immune cells are first, followed by stromal, followed by tumor
    dfClusterT = dfClusterT.reindex(["countA%", "countB%", "countK%",'countC%','countD%','countJ%','countX%','countE%','countN%','countF%','countH%'])
    
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
    layout = {'title':'Average ⍺SMA+ Mesenchymal Cell Neighborhood Cluster Composition','xaxis':{'title':'Cluster'},'yaxis':{'title':'Fraction Present'},'barmode':'stack'}
    fig = {'data':traces,'layout':layout}
    pio.write_image(fig,cwd+'/figures/figure4B.png',format='png',scale=2) #saves figure to 'figures' folder as a png file

    print("Figure 4B saved to 'figures' folder.")



    ##FIGURE 4C - SURVIVAL CURVE, CLUSTER 2

    print('\nFigure 4C')
    
    #calculate abundance of each cluster in each tumor region, save data to csvs
    print('Calculating counts of clusters per tumor region...')
    clusterCountPerROI('dfNeighClusteredN45k6')
    print("Counts calculated. Csvs saved to 'dfCreated' folder.")
    
    #run kaplan function to generate survival curve and return plot to save
    plt = kaplan('dfClustCountsN45k6_avg','clust1PercList') #NOTE: title will say cluster 2, this is to match naming in figure 4B but data originally was called cluster 1
    
    plt.savefig(cwd+'/figures/figure4C.png',format='png') #saves figure to 'figures' folder as a png file
    #plt.show()
    plt.close()
    print("Figure 4C saved to 'figures' folder.")        



    ##FIGURE 4D - SURVIVAL CURVE, CLUSTER 6

    print('\nFigure 4D')
    
    #calculate abundance of each cluster in each tumor region, save data to csvs
    #run kaplan function to generate survival curve and return plot to save
    #this assumes clusterCouontPerROI() function has been run for Figure 4C
    plt = kaplan('dfClustCountsN45k6_avg','clust5PercList') #NOTE: title will say cluster 6, this is to match naming in figure 4B but data originally was called cluster 5
    
    plt.savefig(cwd+'/figures/figure4D.png',format='png') #saves figure to 'figures' folder as a png file
    #plt.show()
    plt.close()
    print("Figure 4D saved to 'figures' folder.")        



    ##FIGURE 4E - SURVIVAL CURVE, CLUSTER 3

    print('\nFigure 4E')
    
    #calculate abundance of each cluster in each tumor region, save data to csvs
    #run kaplan function to generate survival curve and return plot to save
    #this assumes clusterCouontPerROI() function has been run for Figure 4C
    plt = kaplan('dfClustCountsN45k6_avg','clust2PercList') #NOTE: title will say cluster 3, this is to match naming in figure 4B but data originally was called cluster 2
    
    plt.savefig(cwd+'/figures/figure4E.png',format='png') #saves figure to 'figures' folder as a png file
    #plt.show()
    plt.close()   
    print("Figure 4E saved to 'figures' folder.")
    
    
    


if __name__=="__main__":
    fig1()
    fig2()
    fig3()
    fig4()