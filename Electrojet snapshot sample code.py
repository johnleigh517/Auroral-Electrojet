# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 11:13:16 2018

@author: john_
"""




import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import datetime
import scipy.signal
from mpl_toolkits.basemap import Basemap
import math
import matplotlib

days=2  #no of files to read


barmax=800  #where barmax is the max value for the colour bar
B_interval=15  #the gradient of dB/dt in mins


fig= plt.figure(1) #needed to name plots
#creating 4 axes
ax = plt.subplot2grid(shape=(2,6), loc=(0,1), colspan=2)
ax1 = plt.subplot2grid((2,6), loc=(0,3), colspan=2)
ax2 = plt.subplot2grid((2,6), loc=(1,3), colspan=2)

ax3 = plt.subplot2grid((2,6) ,loc= (1,1), colspan=2)

#adding axes for colourbars
cax, _= matplotlib.colorbar.make_axes(ax)
cax1, _= matplotlib.colorbar.make_axes(ax1)
cax3, _= matplotlib.colorbar.make_axes(ax3)


for_save=1000   #used to save image of snapshot, 

for time in range(0,1440*days,5):   #creating a loop of videos
    
    m=Basemap(projection='mill', llcrnrlat=33,llcrnrlon=-15,urcrnrlat=80,urcrnrlon=45, resolution='c',ax=ax)
    
    #setting miller projection to only show Europe on the map
    
    
    m.drawcoastlines()
    #m.bluemarble()   #include if you want colour
       
    
    m.drawparallels(np.arange(10,90,20),labels=[1,1,0,1])
        # draw meridians
    m.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])
    
    

    sites_obs = [ 'ABK','BDV','BFO','BOX','CLF',
             'EBR','ESK','HAD','HLP','HRN','KIV','LER','LON','LVV',
             'LYC','NUR','PAG','SFS','SOD','SPG','SPT','SUA',
             'UPS','VAL','WNG'] #sites to evaluate
    

    
    folder = "/Users/john_/Documents/Auroral_electrojet_storms/Patricks_day_storm/" #folder files are in
    # End INPUTS
    
    
    
    count_0 = 0    #reading parameters from files
    for k1, k2 in enumerate(sites_obs):
        k2_l = k2.lower()
        count = 0
        filenames = sorted(os.listdir(folder))
        filenames = [x for x in filenames if str(k2_l) in x]
    
        for i in filenames:
            fff = folder + i
            f = open(fff, 'r')
            data = f.readlines()
            f.close()
        
            # Read the components of the tensors
            for index, line in enumerate(data):
                if line.startswith('DATE'):
                    val = index
                if str(k2) + 'X' in line:
                    corr = 0
                if str(k2) + 'D' in line:
                    corr = 1
        
            df = pd.read_csv(fff, sep = '\s+', skiprows = int(val)) # Check the rows that need to be skiped
        
            count = count + 1
            print(fff)
        
            if count == 1:
                df_t = df
                print(df_t.shape)
            else:        
                df_t = np.row_stack((df_t, df))
                print(df_t.shape)
        
        
        # Define the datetime variables
        Date = df_t[:,0]
        Hour = df_t[:,1]  
        dates = Date + ' ' + Hour
        dtim = [datetime.datetime.strptime(x, '%Y-%m-%d %H:%M:%S.000') for x in dates]
    
        # Selecte the magnetic data   
        Bx_1 = pd.to_numeric(df_t[:,3])
        By_1 = pd.to_numeric(df_t[:,4])
        Bz = pd.to_numeric(df_t[:,5])
    
        if corr == 0:
            Bx = Bx_1
            By = By_1
            
        if corr == 1:
            Bx = Bx_1 * np.cos(By_1*3.141592/(180*60))
            By = Bx_1 * np.sin(By_1*3.141592/(180*60))
    
        if count_0 == 0:
            B_fields = np.zeros([len(sites_obs),len(Bx),3])
            count_0 = 1
            
        B_fields[k1,:,0] = Bx
        B_fields[k1,:,1] = By
        B_fields[k1,:,2] = Bz
    
        #----------------------------------------------------------------------------
        
        #setting time bars
        min_time=[]
        
        for i in range(1440*days):  #sets up a time array in minutes to plot withB field values
            
            min_time.append(i)
           
        #setting time bars
        
        xb=[min_time[time],min_time[time]]    
    
        
        yb=[-5000,60000]
        
    
        B_field_tot=((B_fields[k1,:,0].T)**2+(B_fields[k1,:,1].T)**2)**(1/2)
        
        
        ax2.plot(min_time,B_field_tot, label = k2)
        ax2.plot(xb,yb,color='black')
        ax2.set_ylim(6000,24000)
        ax2.set_xlabel('Time(minutes)')
        
        ax2.set_ylabel('Bxy [nt]')
        ax2.annotate(str(k2), xy=(1440*days,B_field_tot[1440*days-1]))
    
        
    
    
        #-------------------------------------------------
    
        # Detrend to only ofcus on the variations of the magnetic field    
        B_fields[k1,:,:] = scipy.signal.detrend(B_fields[k1,:,:], axis = 0)
    
    
        B_fieldx=B_fields[k1,:,0].T[time]
        B_fieldy=B_fields[k1,:,1].T[time]
        B_fieldz=B_fields[k1,:,2].T[time]
            
        B_field_total=math.sqrt((B_fieldx**2+B_fieldy**2+B_fieldz**2)) 
            
        dtim_new = dtim[time]
            
        xs=[]
        ys=[]
        s=200  
        
        
        B_field_new=math.sqrt(B_fieldx**2+B_fieldy**2)
        #change B_field here
        dBx=(B_fieldx-B_fields[k1,:,0].T[time-B_interval])/B_interval
        dBy=(B_fieldy-B_fields[k1,:,1].T[time-B_interval])/B_interval
        dBz=(B_fieldz-B_fields[k1,:,2].T[time-B_interval])/B_interval
        
        delta_x=dBx+dBy
        
        if delta_x >=0:
            
            dB=math.sqrt(dBx**2+dBy**2)
        
        if delta_x <0:
            dB=-math.sqrt(dBx**2+dBy**2)
            
      
        
        
        dB1=[]
        
        dB1.append(dB)
        #--------------------------------------------------
        #adding bits for colourbar
        cm = plt.cm.get_cmap('RdBu_r')  #seismic is good colour scale also
        cm_xy=plt.cm.get_cmap('OrRd')
        B=[]
        B.append(B_field_new)
        
        normalize_xy = matplotlib.colors.Normalize(vmin=0, vmax=barmax)
        #setting min & max for colorbar, as min max of B_field
        colors = [cm_xy(normalize_xy(value)) for value in B]
        #setting colours of points to change as B_field_new changes
        B_field_new1=B_fieldz
        B1=[]
        B1.append(B_field_new1)
        
        
        normalize = matplotlib.colors.Normalize(vmin=-barmax, vmax=barmax)
        normalizedB = matplotlib.colors.Normalize(vmin=-15, vmax=15)
        #setting min & max for colorbar, as min max of B_field
        colors1 = [cm(normalize(value)) for value in B1]
        #setting colours of points to change as B_field_new changes
        colors2=[cm(normalizedB(value)) for value in dB1]
        
        #------------------------------------------------------------------------
        
         #now adding arrows to B_h field plot
        x_len=[]
        y_len=[]
        
        
        
        
        if B_fieldx>=0 :
            
            
            
            
            xlen=1/(math.sqrt(1+(B_fieldx/B_fieldy)**2)) #used to keep arrows same length
            
            
        
        if B_fieldx < 0:
            
            
            
            xlen=-1.0/(math.sqrt(1+(B_fieldx/B_fieldy)**2))
            
            
        if B_fieldy>=0:
            
            
            
            ylen=1/(math.sqrt(1+(B_fieldy/B_fieldx)**2))    
            
        
        if B_fieldy < 0:
            
            
            ylen=-1.0/(math.sqrt(1+(B_fieldy/B_fieldx)**2))
            
        
        dx=400000*xlen   #sets length of arrow
        dy=400000*ylen
        


        #parametersof arrows
        linewidth=2.5
        head_width=150000
        head_length=100000
         
    #-----------------------------------------------------
            
    
        if k2 == 'ABK'or'BDV'or'BEL'or'BFO'or'BOX'or'CLF'or'DOU'or'EBR'or'ESK'or'FUR'or'GCK'or'HAD'or'HLP'or'HRB'or'HRN'or'KIV'or'LER'or'LON'or'LVV'or'LYC'or'MAB'or'NCK'or'NGK'or'NUR'or'PAG'or'SFS'or'SOD'or'SPG'or'SPT'or'SUA'or'THY'or'UPS'or'VAL'or'WIC'or'WNG':
          #reads whether a file is present and divides each file type
            
            
            
            #for index, line in enumerate(data):     
            xpt=data[5]   #reading longitude
            x=[]
            x.append(xpt)   #changing line to a list
        
            xpt_t=xpt.split()   #seperating list
            xpt_t1=xpt_t[2]        #Taking numerical value from line
            
            xpt_0=float(xpt_t1)   #converting numerical value from string to float
            
            if xpt_0 > 180:
                xpt_0=xpt_0-360  #needed to correct for points in western europe
            else:
                xpt_0=xpt_0
            
        
            ypt=data[4]    #reading latitude from files
            y=[]
            y.append(ypt)
        
            ypt_t=ypt.split()
            ypt_t1=ypt_t[2]
            ypt_0=float(ypt_t1)
            
            xpt_prop,ypt_prop=m(xpt_0,ypt_0)
            
    
            xs.append(xpt_prop)
            ys.append(ypt_prop)
            ax.scatter(xs,ys,s=s,color=colors,edgecolor='black')
            ax.arrow(xpt_prop,ypt_prop,dx,dy,fc="k", ec="k", linewidth = linewidth, head_width=head_width,head_length=head_length)
            ax1.scatter(xs,ys,s=s,color=colors1,edgecolor='black')
            n=[k2] #needing to define name of magnetometer
            ax3.scatter(xs,ys,s=s,color=colors2,edgecolor='black')
            for i, txt in enumerate(n):
                #ax.annotate(txt, (xs[i], ys[i]))
                ax1.annotate(txt, (xs[i], ys[i]))
    
    
    
    
    #creating second & third basemap for Bz and dB/dt
    
        m1=Basemap(projection='mill', llcrnrlat=33,llcrnrlon=-15,urcrnrlat=80,urcrnrlon=45, resolution='c',ax=ax1)
        
        #setting miller projection to only show Europe on the map
        
        
        m1.drawcoastlines()
    
        
        m1.drawparallels(np.arange(10,90,20),labels=[1,1,0,1])
            # draw meridians
        m1.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])
        #m1.etopo()
        m2=Basemap(projection='mill', llcrnrlat=33,llcrnrlon=-15,urcrnrlat=80,urcrnrlon=45, resolution='c',ax=ax3)
        
        #setting miller projection to only show Europe on the map
        
        
        m2.drawcoastlines()
    
        
        m2.drawparallels(np.arange(10,90,20),labels=[1,1,0,1])
            # draw meridians
        m2.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])
       
     
    
    ax.set_title( str(dtim_new) +'\nAbsolute Magnitude of Horizontal Field')
    ax1.set_title( str(dtim_new) +'\nMagnitude of Z-Field')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')    
    ax3.set_title('dBx/dt(nT/min)')   
    #ax2.set_title('Horizontal B field Strength')
    
    #adding colourbars
    cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cm_xy, norm=normalize_xy,label= 'B Field Variation[nT]')
    
    
    cbar = matplotlib.colorbar.ColorbarBase(cax1, cmap=cm, norm=normalize,label= 'B Field Variation[nT]')
    cbar= matplotlib.colorbar.ColorbarBase(cax3, cmap=cm, norm=normalizedB)
    
    plt.pause(0.5) 
    
    #---------------------------------------------------
    #code used to save each image for each loop
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.show()
    plt.savefig('/Users/john_/Documents/Auroral Electrojet Images/Patricks_day_storm/'+str(for_save) +'.png')
    for_save=for_save+5
    
    
    #---------------------------------------------
    
    #clearing axes to allow loop to run again, without overplotting
    ax.clear()
    ax1.clear()
    ax2.clear()

    ax3.clear()
    