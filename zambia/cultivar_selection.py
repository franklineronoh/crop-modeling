# -*- coding: utf-8 -*-
"""
Created on Mon May  1 11:51:33 2023

author: Sara Miller
calculate district level stats for RHEAS vs preharvest survey data for Zambia
select best cultivars for each district
"""

import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import sys
import glob
import re

# set working directory
wdir = r'D:\RONO\rheas\maize_crop_modeling\crop-modeling\zambia\data\cultivartests\\'

###############################################################################
r'''
# first, loop through all CSV files containing the different cultivars and years
# then append all data to one CSV

# get list of csvs containing cultivar test results
# these files should be in a sub-folder in your working directory called 'csvs'
flist = glob.glob(wdir+'cultivars\cultivar*.csv')
df = pd.DataFrame()
for f in flist:
    df1 = pd.read_csv(f)
    # get cultivar number and year from file name, add to csv
    cultivar = re.search(r'cultivars\\cultivar(.*)_',f).group(1)
    df1['cultivar'] = cultivar
    year = re.search(r'cultivars\\cultivar{0}_(.*).csv'.format(cultivar),f).group(1)
    df1['year'] = int(year)
    # make sure that only yields from the selected year are added
    df1['planting'] = pd.to_datetime(df1['planting'])
    df1 = df1.set_index('planting')
    df1['PlantingYear'] = df1.index.year.values
    df1['PlantingMonth'] = df1.index.month.values
    df1 = df1.reset_index()
    df1 = df1.loc[((df1['PlantingYear']==int(year))&(df1['PlantingMonth']>8))]
    # append to one single dataframe and save to a CSV
    df = pd.concat([df,df1],ignore_index=True)

    
df.to_csv(wdir+'allcultivars.csv', index=False)'''

###############################################################################
r'''
# read in csvs with dssat outputs, survey data
pre = pd.read_csv(r'D:\RONO\rheas\maize_crop_modeling\crop-modeling\zambia\data\preharvest\zambia_preharvest_district_data_clean.csv')
df = pd.read_csv(wdir+'allcultivars.csv')


# optionally: remove anomalous yield values from the preharvest dataset before comparing to DSSAT
#pre = pre.loc[pre['flag']=!2]

df = df.rename(columns={'cname':'district'})
pre = pre.loc[((pre['year']>=2011)&(pre['year']<2022))]

# convert columns to lower and cname names
df.columns=df.columns.str.lower()
df.district=df.district.str.lower()
df.district=df.district.str.replace('-', '', regex=True)
df.district=df.district.str.replace(' ', '', regex=True)



# df = df.reset_index()
pre = pre.reset_index()


stats = pd.DataFrame()

for c in list(df.cultivar.unique()):
    for p in list(pre.district.unique()):
        pre1 = pre.loc[pre['district']==p].copy()
        df2 = df.loc[((df['district']==p)&(df['cultivar']==c))].copy()

        # get average season length in days for each cultivar
        try:
            df2['harvest'] = pd.to_datetime(df2['harvest'])
            df2['planting'] = pd.to_datetime(df2['planting'])
            
            df2['slength'] = (df2['harvest']-df2['planting']).dt.days
            slength = df2['slength'].mean()
        except:
            slength = np.NaN

        
        # pre1 = pre1.set_index(['district','year'])
        df2 = df2.groupby(['district','year']).mean()
        pre1 = pre1.groupby(['district','year']).mean()

        # pre1 = pre1.set_index(['district','year'])
        # df2 = df2.set_index(['district','year'])


        
        # combine the two dataframes: note that here the yield is stored in column 'gwad' for dssat
        # and column 'yield' for the survey data
        # note that if there are years that one dataset has data for that the other doesn't
        # the years that are not present will be filled with NaNs
        df1 = pd.concat([pre1,df2],axis=1)

        # calculate statistics one province at a time and append to a new dataframe
        # correlation
        corr = df1['gwad'].corr(df1['yield'])
        # root mean square error
        rmse = np.sqrt(((df1['gwad']-df1['yield'])**2).mean())
        # percentage rmse
        perrmse = np.sqrt(((df1['gwad']-df1['yield'])**2).mean())/(df1['yield'].mean())*100
        # unbiased rmse
        ubrmse = np.sqrt((((df1['gwad']-df1['gwad'].mean())-(df1['yield']-df1['yield'].mean()))**2).mean())
        # bias
        bias = (df1['gwad']-df1['yield']).mean()
        stats1 = pd.DataFrame({'district':[p],'cultivar':[c], 'rmse':[rmse], '%rmse':[perrmse],'unbiased rmse':[ubrmse], 'bias':[bias], 'correlation':[corr], 'season length':[slength]})
        stats = pd.concat([stats,stats1],ignore_index=True)
        # stats = stats.append({'district':p,'cultivar':c, 'rmse':rmse, '%rmse':perrmse,'unbiased rmse':ubrmse, 'bias':bias, 'correlation':corr, 'season length':slength}, ignore_index=True)

# save statistics as a csv
stats.to_csv(wdir+'districteval.csv')

###############################################################################
# select best cultivars for each district given the statistics
culs = pd.read_csv(wdir+'districteval.csv')
selected = pd.DataFrame()

for p in list(culs.district.unique()):

        df = culs.loc[culs['district']==p]
        # remove cultivars that have a season length shorter than expected
        earliestharvest = 100 # shortest expected season length in days
        df = df.loc[df['season length']>=earliestharvest]
        # select cultivars with correlation above 0.8 and the lowest ubrmse
        df1 = df.loc[df['correlation']>=0.8]
        if df1.empty:
            # if no cultivars have correlation above 0.8 select the cultivar with highest correlation
            df1 = df.loc[df['correlation'] == df['correlation'].max()]
        else: 
            df1 = df1.loc[df1['unbiased rmse'] == df1['unbiased rmse'].min()]
        selected = pd.concat([selected,df1])

# save stats and cultivar number of selected cultivars to a csv
selected.to_csv(wdir+'cultivartests_selected.csv')

# add the dssat outputs from only the selected cultivars to one csv
yields = pd.read_csv(wdir+'allcultivars.csv')
selected = pd.read_csv(wdir+'cultivartests_selected.csv')
yields2 = pd.DataFrame()

# convert columns to lower and cname names
yields.columns=yields.columns.str.lower()
yields.cname=yields.cname.str.lower()
yields.cname=yields.cname.str.replace('-', '', regex=True)
yields.cname=yields.cname.str.replace(' ', '', regex=True)

for c in list(selected.district.unique()):
    sn = selected.loc[(selected['district']==c)]
    df = selected.loc[((selected['district']==c))]
    y1 = yields.loc[((yields['cname']==c))]
    cul = df['cultivar'].item()
    y2 = y1.loc[y1['cultivar']==cul]
    yields2 = pd.concat([yields2,y2])

yields2.to_csv(wdir+'cultivarselectedyields.csv')'''

###############################################################################

# plot statistics of the selected cultivars
# read in shapefiles and statistics of selected cultivars
zmb2 = gpd.read_file(r'D:\RONO\rheas\maize_crop_modeling\crop-modeling\zambia\data\shp\zmb_province_cleaned.shp')
zmb = gpd.read_file(r'D:\RONO\rheas\maize_crop_modeling\crop-modeling\zambia\data\shp\zmb_district_cleaned.shp')
stats = pd.read_csv(wdir+'cultivartests_selected.csv')
stats.rename(columns={'district':'district_b'}, inplace=True)

# rename province column to match with shapefile province name column for plotting
# note: you may have to edit some names if there are inconsistencies in spelling or capitalization
# stats = stats.rename(columns={'district':'NAME1'})
zmb = zmb.merge(stats, on='district_b', how='left')

# plot correlation, bias, and unbiased RMSE at province level
fig,ax = plt.subplots(1,3,figsize=(8,3))
zmb.plot(column='correlation',cmap='BrBG',edgecolor='black',legend=True,ax=ax[0])
zmb.plot(column='bias',cmap='Reds',edgecolor='black',legend=True,ax=ax[1])
zmb.plot(column='unbiased rmse',cmap='YlOrBr',edgecolor='black',legend=True,ax=ax[2])
ax[0].set(title='Correlation')
ax[1].set(title='Bias')
ax[2].set(title='Unbiased RMSE')
plt.tight_layout()
plt.show()

