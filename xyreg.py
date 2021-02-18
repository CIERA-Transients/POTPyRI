#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 21:47:32 2021

@author: dvelasco
"""
    
#To make reg file from test.cat from sextraxt

def xyreg(xpos, ypos, fname, color):
    f= open("/Users/dvelasco/Imaging_pipelines/"+fname,"w+")  
    line1 = '# Region file format: DS9 version 4.1'
    line2 = 'global color='+color+' dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
    line3= 'PHYSICAL'
    #print(len(stellae_ra))
    f.write(line1+"\n")
    f.write(line2+"\n") 
    f.write(line3+"\n") 
    for i in range(len(xpos)):
        coords = 'point('+str(xpos[i])+','+str(ypos[i])+') # point=circle'
        f.write(coords+"\n") 
    f.close()