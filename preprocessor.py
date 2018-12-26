#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 10:44:56 2018

@author: konrad
"""
from __future__ import division
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from math import sqrt,pi,exp, log
import os
import re
import urllib

def download_links(token, dest_dir):
    url = "http://refdata:research@62.218.164.232/lva2019/" + token
    url_list = []
    
    append_list = ["_annotation.csv", "_annotation_user.csv", "_location.csv", "_sensor.csv"]
    
    current_dir = os.getcwd()
    
    os.chdir(dest_dir)
    
    try:
        url_file = urllib.urlopen(url)
        url_text = url_file.read()
        
        matches = re.findall(token+'_\d*-\d*/', url_text)
        if matches:
            #print "Found Match!"
            for match in matches:
                #print match
                url_list.append(url+"/"+match+match)
    
    except IOError:
        print "Problem: ", IOError
        
    
    #Get all the download links
    for link in url_list:       
        #Download the trips    
        for append_word in append_list:
            #Generate the links for downloading
            print link[0:-1] + append_word
            urllib.urlretrieve(link[0:-1] + append_word, link[-32:-1]+append_word)
            
    os.chdir(current_dir)

def PAA(time,array):
    """
    Does a piecewise aggregate function on the array
    """
    Sum = 0
    num_of_items = 0
    current_time = time[0]
    new_time = []
    new_array = []
    
    for i in range(len(time)):
        #aggregate every 50 miliseconds
        if(time[i]-current_time > 50):
            new_array.append(Sum/num_of_items)
            new_time.append(current_time + (time[i]-current_time)/2)
            
            current_time = time[i]
            Sum = 0
            num_of_items = 0
            
        #do the average
        Sum += array[i]
        num_of_items += 1
        
    return np.array(new_time),np.array(new_array)
    

#function to get the array of acceleration
def interpolate_accelearation(dframe):
    tp = np.array(dframe['time'])
    x = np.array(dframe['x'])
    y = np.array(dframe['y'])
    z = np.array(dframe['z'])
    
    #interpolation with linear splines
    inter_x = interp1d(tp, x, kind='slinear')
    inter_y = interp1d(tp, y, kind='slinear')
    inter_z = interp1d(tp, z, kind='slinear')
    
    time_range = np.arange(tp[0],tp[-1])
    
    #interpolation for data
    ax = inter_x(time_range)
    ay = inter_y(time_range)
    az = inter_z(time_range)
    
    #Piecewise Aggregate Function
    _,ax = PAA(time_range,ax)
    _,ay = PAA(time_range,ay)
    time_range,az = PAA(time_range,az)
    
    time_range -= tp[0]
    time_range /= 1000
    
    #calculate total acceleration
    total = []
    for x,y,z in zip(ax,ay,az):
        total.append(sqrt(x**2+y**2*z**2))
        
    total = np.array(total)
    
    names = ["time","x-acceleration","y-acceleration","z-acceleration","total-acceleration"]
    indices = range(len(time_range))
    
    df = pd.DataFrame(np.array([time_range,ax,ay,az,total]).transpose(), index = indices, columns = names)
    
    df['mode'] = dframe['mode']
    df['notes'] = dframe['notes']

    return df
        
def get_filelist(dir, key):
    """
    returns a list with the corresponding csv names of type *key.csv
    """
    filenames = os.listdir(dir)
    
    l = []
    
    for name in filenames:
        match = re.search(key, name)
        if match:
            l.append(dir+name)
            
    return l

def fast_fourier_transform(df, key, n_recon = 30):
    """
    Get and visualize Fourier Transform
    key in {"x-acceleration","y-acceleration","z-acceleration","total-acceleration"}
    n_recon number of Spectrals/Frequencies used for reconstruction of signal
    """
    
    Spektrum = np.fft.fft(df[key])
    freq = np.fft.fftfreq(len(df['time']))
    
    plt.figure(1)
    plt.subplot(211)
    plt.plot(freq, Spektrum.real)
    #plt.xscale('log')
    plt.yscale('log')
    
    plt.subplot(212)
    plt.plot(freq, Spektrum.imag, 'r')
    #plt.xscale('log')
    plt.yscale('log')
    
    #amplying inverse of fourier
    #TODO: Scaling of Signal after reconstruction
    inv = np.fft.ifft(Spektrum[0:n_recon])
    plt.figure(2)
    plt.subplot(211)
    plt.plot(inv)
    
def split_in_segments(df):
    """
        Splits the data frame in 30 seconds groups
        returns a list that contains the DataFrames of 30 seconds
    """
    #list that holds 4 lists (x, y, z, n2) of 30s blocks
    lst = []
    
    idx = 0
    time = []
    xlst = []
    ylst = []
    zlst = []
    n2lst = []
    
    for i in range(len(df['time'])):
        #check if a new 30s block must be created
        if i - idx >= 600:
            idx = i
            
            seg = np.array([time, xlst, ylst, zlst, n2lst]).transpose()
            dframe = pd.DataFrame(seg, columns = df.columns[0:-2])
            dframe['mode'] = df['mode']
            dframe['notes'] = df['notes']
            lst.append(dframe)
            
            time = []
            xlst = []
            ylst = []
            zlst = []
            n2lst = []
            
        #Current row of dataframe
        row = df.iloc[i]
        
        time.append(row[0])
        xlst.append(row[1])
        ylst.append(row[2])
        zlst.append(row[3])
        n2lst.append(row[4])
        
    return lst
        
#Token
mytoken = "357810086663607"

#Directory path to stored trips relative to working directory!!!
directory_name = "/home/konrad/Documents/UniWien/WS18/Data Mining/Test/"
word = "sensor.csv"

download_links(mytoken, directory_name)

#get List of all relevant filenames
sensor_files = get_filelist(directory_name,"sensor.csv")
annotations_names = get_filelist(directory_name,"annotation.csv")

index = 3

sensor = sensor_files[index]
annotation = annotations_names[index]

#read in the dataframe
sensordata = pd.read_csv(sensor)
annotationdata = pd.read_csv(annotation)

sensordata['mode'] = annotationdata['mode'][0]
sensordata['notes'] = annotationdata['notes'][0]

#describe the dataframe
sensordata.describe()

#erase duplictes 
sensordata = sensordata[sensordata['sensor']=="acceleration"]
sensordata.drop_duplicates("time", inplace=True)

acceleration = interpolate_accelearation(sensordata)

acceleration.plot(x="time",y="total-acceleration")
fast_fourier_transform(acceleration, 'total-acceleration',50)
