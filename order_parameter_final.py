"""
Created on Sun Aug 26 23:12:43 2018
@author: harris-a
"""
#Purpose: To write a program that calculates the average number of neighbors to 
#           a group of particles when given multiple text files with coordinates.
#Process: Feeds in txt files, converts to a list of dataframes, calculates the 
#           distances between each particle, counts the number of neighbors for
#           each particle in a frame and then creates a data frame with the 
#           number of neighbors per particle and a list with the average per frame. 
#Output: A list of the average number of particle neighbors (first order, second order, 
#        etc.) per frame, list of the calculated number of neighbors (first order,
#        second order, etc.) per particle in a frame, and histograms of the
#        the number of first order neighbors per particle in a frame. 
#Developed by: Aaron Harris and Alicia Altemose 8/27/2018
 

###############################################################################
#Variable Editor- Change run parameters
###############################################################################
#location of the .txt files
path = 'C:/Users/User/Documents/Coord-files/' 

#Minimun and Maximum Distance Range for Neighbors 
min_ = 0                                     
max_ = 6                                     
inclusive = False 


zero= False                           

###############################################################################
#Neighbors Calculator Program
###############################################################################
#Import Libraries
import pandas as pd
import numpy as np
from glob import glob

#Create the function for calculating the distance and determining neighbors. 
def calc_neighbors(cord_df,min_range=0, max_range=6, inclu=False):
    '''
    This function is created to feed in a dataframe of x and y coordinates along with
    the minimum and maximum distance range for a neighbor. The function will calculate 
    the distances, determine if a record is a neighbor or not for the set.
    This neighbor count per particle is written to a dataframe for further analysis. 
    
    Argument Definitions:
        cord_df= A dataframe with no headers and the variables of x and y in that order (type=pd.DataFrame)
        min_range= The lower limit of what is defined as a neighbor. Default=0 (type=float)
        max_range= The upper limit of what is defined as a neighbor. Default=6 (type=float)
        inclu= Argument to specify if neighbor range is inclusive or exclusive (type=bolean)
    '''
    #set up initial values to loop through in tuples and list. 
    cord_tups = [tuple(x) for x in cord_df.values]  
    x_ls = list(cord_df.iloc[:,0]) 
    y_ls = list(cord_df.iloc[:,1])  
    dist_df = pd.DataFrame()                 
    #For loop that measures the distance between each particle and writes it to 
       #a data frame (matrix) of each particle against another. 
    for x,y in cord_tups:
        dist_ls = []
        for i in range(0,len(cord_df),1):
            x_dist = x-x_ls[i]
            y_dist = y-y_ls[i]
            dist_ls.append((x_dist**2 + y_dist**2)**0.5)
        dist_df = dist_df.append(pd.Series(dist_ls), ignore_index=True)
    #counts the number of neighbors 
    return ([pd.Series.between( dist_df.iloc[:,n], min_range,
                               max_range, inclu).sum() for n in dist_df.columns],dist_df)


def calc_next_neighbors(distance_df, nei_ls, n_min= 6, order=1, partic=6):
    '''
    Function to count neighbors greater than the first order. This function can be used
    for 2 to N orders of neighbors. This function returns a list of the count of neighbors
    at the given order. 
    
    Argument Definitions:
        distance_df = A pandas data frame with distances between each particle. 
        nei_ls = A list of the number of neighbors in the order below the order being calculated. 
        n_min = The size (diameter) of the particles being analyzed. 
        order = One less than the order of neighbors being calculated (1 for second order). 
        partic = The minimum number of neighbors for the order to be calculated.
    '''
    sec_nei_ls =[]
    for y,x in zip(range(len(nei_ls)),nei_ls):
        if x >= partic*order:
            sec_nei = distance_df[y][(distance_df[y] >= n_min*order) & 
                                 (distance_df[y] <= n_min*(order+1))].count()
            if sec_nei > partic*(order+1):
                sec_nei = partic*(order+1)
            else:
                sec_nei = sec_nei
        else:
            sec_nei = 0
        sec_nei_ls.append(sec_nei)
    return sec_nei_ls
                  
    
#Create a function to plot the distributions
def plot_neighbor_dis(n_ls, include_zero=True):
    '''
    Function to create a list of histograms to display the number of neighbors in
    each frame. 
    
    Argument Defintions:
        n_ls = A list of the number of neighbors in a given frame. 
        include_zero = Boolean identifier for if zero neighbors should be included in the graph. 
    '''
    if include_zero == True:
        n_ls = pd.Series(n_ls).value_counts()
        n_plot=n_ls.plot.bar(width=.85,figsize=(5,5),
                             title='Distribution of Neighbors')
    else:
        n_ls = pd.Series(n_ls)
        n_ls = n_ls[n_ls !=0]
        n_ls = pd.Series(n_ls).value_counts()
        n_plot= n_ls.plot.bar(width=.85,figsize=(5,5),
                              title='Distribution of Non-Zero Neighbors')
    return n_plot


# Import the text files from the folder into a list of data frames to 
#  iterate through. 
Frames = glob( path + "*.txt")
df_ls= []
for f in Frames:
    df = pd.read_csv(f,index_col=None, sep="\t", header=None)
    df_ls.append(df)

#Create the matrix of the number of neighbors per particle and the distances between
# each particle, along with the list of average numbers of neighbors per frame. 
neighbors_ls= []
avg_ls=[]
avg_no0_ls=[]
dis_matrix_ls=[]
for x in df_ls:
    cn, dis= calc_neighbors(x,min_, max_, inclusive)
    avg= sum(cn)/len(cn)
    avg_no0= np.mean(list(filter(lambda x: x != 0, cn)))
    dis_matrix_ls.append(dis)
    neighbors_ls.append(cn)
    avg_ls.append(avg)
    avg_no0_ls.append(avg_no0)
    
    
#Calculate the number of second order neighbors, along with the average of those per frame. 
sec_neighbors_ls= []
sec_avg_ls=[]
sec_avg_no0_ls = []
for y,x in zip(range(len(neighbors_ls)),neighbors_ls):
    cn2= calc_next_neighbors(dis_matrix_ls[y], x)
    avg2= np.mean(list(filter(lambda x: x >= 0, cn2)))
    avg2_no0= np.mean(list(filter(lambda x: x != 0, cn2)))
    sec_neighbors_ls.append(cn2)
    sec_avg_ls.append(avg2)
    sec_avg_no0_ls.append(avg2_no0)
    
#Creates a list of histograms (one for each frame) of the number of neighbors per particle.     
dis_graphs= []
for i in range(len(neighbors_ls)):
    dis_graphs.append([plot_neighbor_dis(neighbors_ls[i],zero)])
    i = i+1