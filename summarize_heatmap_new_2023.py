#!/usr/bin/python
import os,sys,glob,re
import numpy as np
import scipy
from scipy import stats
import datetime
import time
from datetime import timedelta
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from matplotlib import colors as c
#from matplotlib  import cm
from scipy.stats.kde import gaussian_kde
from numpy import linspace
from scipy.stats import kruskal
#from scipy.stats import nanmean
#from scipy.stats import nanmedian
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mstats


#_responselatency_: Response latency (ms) 
#_responsefrequency_: Response frequency

#_responsecumulativedistance_: Response cumulative distance (pixels)
#_responsecumulativedpix_: Response cumulative dpix (pixels)
#_responsedisplacement_: Response displacement (pixels)
#_responsepeakdpix_: Peak dpix (pixels)
#_responsepeakspeed_: Peak speed (pixels / ms)
#_responsepolygonarea_: Response movement area (pixels)
#_responsespeed_: Response speed (pixels / ms)
#_responseaveinitangvel_: Response ave angular velocity (degrees / ms)
#_responsetime_: Response time (ms)
#_responsevelocity_: Response velocity (pixels / ms)

#habituation_day5dpfhab1pre=1_15:35:00-1_15:54:00
#habituation_day5dpfhab1=1_16:03:00-1_16:05:00
#habituation_day5dpfhab1post=1_16:09:00-1_16:28:00
#habituation_day5dpfhab2wdf=1_16:37:00-1_16:39:00
#habituation_day5dpfhab2post=1_16:43:00-1_17:03:00
#habituation_day5dpfhab3wdf=1_17:11:00-1_17:13:00
#habituation_day5dpfhab3post=1_18:07:00-1_18:26:00
#habituation_day5dpfhab3postshort=1_18:28:00-1_18:34:00
#darkflash_day6dpfdfall=2_9:59:00-2_15:01:00
#darkflash_day6dpfdf1=2_9:59:00-2_11:00:00
#darkflash_day6dpfdf1a=2_9:59:00-2_10:09:30
#darkflash_day6dpfdf1b=2_10:49:30-2_11:00:00
#darkflash_day6dpfdf2=2_11:59:00-2_12:59:30
#darkflash_day6dpfdf2a=2_11:59:00-2_12:09:30
#darkflash_day6dpfdf2b=2_12:49:30-2_12:59:30
#darkflash_day6dpfdf3=2_13:59:00-2_15:01:00
#darkflash_day6dpfdf3a=2_13:59:00-2_14:09:30
#darkflash_day6dpfdf3b=2_14:50:30-2_15:01:00

nonstimcombos = {"Seizure-like movement [day & night]":["_boutseizurecount_"],"Frequency of movement [day]":  ["_numberofboutsSLEEP_","_wakingactive_", "_active_","_numberofbouts_", "interbouttime_"], "Location in well [day]": ["_interboutcenterfraction_","_boutcenterfraction_", "_interboutrhofraction_","_boutrhofraction_"], "Magnitude of movement [day]": ["_boutaveangvelocity_","_boutcumulativemovement_","_boutdisplacement_","_boutpeakangvelocity_", "_boutrevolutions_","_boutspeed_", "_bouttime_", "_boutvelocity_"],"Frequency of movement [night]":  ["_numberofboutsSLEEP_","_wakingactive_", "_active_","_numberofbouts_", "interbouttime_"], "Location in well [night]": ["_interboutcenterfraction_","_boutcenterfraction_", "_interboutrhofraction_","_boutrhofraction_"], "Magnitude of movement [night]": ["_boutaveangvelocity_","_boutcumulativemovement_","_boutdisplacement_","_boutpeakangvelocity_", "_boutrevolutions_","_boutspeed_", "_bouttime_", "_boutvelocity_"]}

skiplist = ["_responsefulldata_", "_responsefulldpixdata_", "_responseinitdirectionheadingangle_", "_responsesumabsheadingangle_", "_responsesumheadingangle_"]
#typecombos = [["Night tap habituation", "Day tap habituation 1", "Day tap habituation 2", "Day tap habituation 3"], ["Day light flash", "Night light flash"],["Night early prepulse tap", "Day early prepulse tap"], ["Night all prepulse tap", "Day all prepulse tap"], ["Day all strong tap", "Night all strong tap"], ["Day early strong tap","Night early strong tap"],["Night early weak tap", "Day early weak tap"], ["Day all weak tap", "Night all weak tap"], ["Dark flash block 3 start","Dark flash block 3 end","Dark flash block 4 start","Dark flash block 4 end","Dark flash block 1 start","Dark flash block 1 end","Dark flash block 2 start","Dark flash block 2 end"]]
#ribgraph_mean_ppi_day6dpfppinight_responsevelocity_1_a001f1000d5pD300a1f1000d5p_a001%4%12.png
stimcombos = {
	"Light flash": [("_day5dpflfday_", "_v1000D1000v200_"),("_day6dpflfnight_","_v1000D1000v200_"),("_day6dpfmslf_", "_v1000D1000v200_")],
	"Tap habituation": [("_day5dpfhab1_","_a1f1000d5p_"), ("_day5dpfhab2wdf_","_a1f1000d5p_"), ("_day5dpfhab3wdf_","_a1f1000d5p_")],
	"Tap post habituation": [("_day5dpfhab3post_","_a1f1000d5p_"), ("_day5dpfhab2post_","_a1f1000d5p_"), ("_day5dpfhab3post_","_a1f1000d5p_")],
	"Prepulse tap": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a1%89%97"),("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a1%89%97")],
	"Weak tap": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day5dpfppi_","_a001f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5p_a001%4%12"),("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day6dpfppinight_","_a001f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5p_a001%4%12")],
	"Strong tap": [("_day5dpfppi_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a0f1400d5pD300a1f1400d5p_a1%89%97"),("_day5dpfhab1pre_","_a1f1000d5p_a1%4%12"),("_day6dpfppinight_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a0f1400d5pD300a1f1400d5p_a1%89%97")],
	"Dark flash": [("_day6dpfdf","_v0D1000v200_v0%18%266")]
	}

stimcombosnorm = {
	"Light flash frequency": [("_day5dpflfday_", "_v1000D1000v200_"),("_day6dpflfnight_","_v1000D1000v200_"),("_day5dpflfday_", "_v1000_"),("_day6dpflfnight_","_v1000_"),("_day6dpfmslf_", "_v1000D1000v200_")],
	"Light flash magnitude": [("_day5dpflfday_", "_v1000D1000v200_"),("_day6dpflfnight_","_v1000D1000v200_"),("_day5dpflfday_", "_v1000_"),("_day6dpflfnight_","_v1000_"),("_day6dpfmslf_", "_v1000D1000v200_")],
	"Light flash latency": [("_day5dpflfday_", "_v1000D1000v200_"),("_day6dpflfnight_","_v1000D1000v200_"),("_day5dpflfday_", "_v1000_"),("_day6dpflfnight_","_v1000_"),("_day6dpfmslf_", "_v1000D1000v200_")],
	"Tap habituation frequency": [("_day5dpfhab1_","_a1f1000d5p_"), ("_day5dpfhab2wdf_","_a1f1000d5p_"), ("_day5dpfhab3wdf_","_a1f1000d5p_")],
	"Tap habituation magnitude": [("_day5dpfhab1_","_a1f1000d5p_"), ("_day5dpfhab2wdf_","_a1f1000d5p_"), ("_day5dpfhab3wdf_","_a1f1000d5p_")],
	"Tap habituation latency": [("_day5dpfhab1_","_a1f1000d5p_"), ("_day5dpfhab2wdf_","_a1f1000d5p_"), ("_day5dpfhab3wdf_","_a1f1000d5p_")],
	"Tap post habituation frequency": [("_day5dpfhab3post_","_a1f1000d5p_"), ("_day5dpfhab2post_","_a1f1000d5p_"), ("_day5dpfhab3post_","_a1f1000d5p_")],
	"Tap post habituation magnitude": [("_day5dpfhab3post_","_a1f1000d5p_"), ("_day5dpfhab2post_","_a1f1000d5p_"), ("_day5dpfhab3post_","_a1f1000d5p_")],
	"Tap post habituation latency": [("_day5dpfhab3post_","_a1f1000d5p_"), ("_day5dpfhab2post_","_a1f1000d5p_"), ("_day5dpfhab3post_","_a1f1000d5p_")],
	"Prepulse tap frequency": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a1%89%97"),("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a1%89%97")],
	"Prepulse tap magnitude": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a1%89%97"),("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a1%89%97")],
	"Prepulse tap latency": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a1%89%97"),("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a1%89%97")],
	"Weak tap frequency": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day5dpfppi_","_a001f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5p_a001%4%12"),("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day6dpfppinight_","_a001f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5p_a001%4%12")],
	"Weak tap magnitude": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day5dpfppi_","_a001f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5p_a001%4%12"),("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day6dpfppinight_","_a001f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5p_a001%4%12")],
	"Weak tap latency": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day5dpfppi_","_a001f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5p_a001%4%12"),("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day6dpfppinight_","_a001f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5p_a001%4%12")],
	"Strong tap frequency": [("_day5dpfppi_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a0f1400d5pD300a1f1400d5p_a1%89%97"),("_day5dpfhab1pre_","_a1f1000d5p_a1%4%12"),("_day6dpfppinight_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a0f1400d5pD300a1f1400d5p_a1%89%97")],
	"Strong tap magnitude": [("_day5dpfppi_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a0f1400d5pD300a1f1400d5p_a1%89%97"),("_day5dpfhab1pre_","_a1f1000d5p_a1%4%12"),("_day6dpfppinight_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a0f1400d5pD300a1f1400d5p_a1%89%97")],
	"Strong tap latency": [("_day5dpfppi_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a0f1400d5pD300a1f1400d5p_a1%89%97"),("_day5dpfhab1pre_","_a1f1000d5p_a1%4%12"),("_day6dpfppinight_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a0f1400d5pD300a1f1400d5p_a1%89%97")],
	"DF end blocks  frequency": [("_day6dpfdf1b_","_v0D1000v200_v0%18%266"), ("_day6dpfdf2b_","_v0D1000v200_v0%18%266"), ("_day6dpfdf3b_","_v0D1000v200_v0%18%266")],
	"DF end blocks  magnitude": [("_day6dpfdf1b_","_v0D1000v200_v0%18%266"), ("_day6dpfdf2b_","_v0D1000v200_v0%18%266"), ("_day6dpfdf3b_","_v0D1000v200_v0%18%266")],
	"DF end blocks  latency": [("_day6dpfdf1b_","_v0D1000v200_v0%18%266"), ("_day6dpfdf2b_","_v0D1000v200_v0%18%266"), ("_day6dpfdf3b_","_v0D1000v200_v0%18%266")],
	"DF begin blocks  frequency": [("_day6dpfdf1a_","_v0D1000v200_v0%18%266"), ("_day6dpfdf2a_","_v0D1000v200_v0%18%266"), ("_day6dpfdf3a_","_v0D1000v200_v0%18%266")],
	"DF begin blocks  magnitude": [("_day6dpfdf1a_","_v0D1000v200_v0%18%266"), ("_day6dpfdf2a_","_v0D1000v200_v0%18%266"), ("_day6dpfdf3a_","_v0D1000v200_v0%18%266")],
	"DF begin blocks  latency": [("_day6dpfdf1a_","_v0D1000v200_v0%18%266"), ("_day6dpfdf2a_","_v0D1000v200_v0%18%266"), ("_day6dpfdf3a_","_v0D1000v200_v0%18%266")],
	"Dark flash all frequency": [("_day6dpfdfall_","_v0D1000v200_v0%18%266"),("_day6dpfdf1_","_v0D1000v200_v0%18%266"),("_day6dpfdf2_","_v0D1000v200_v0%18%266"),("_day6dpfdf3_","_v0D1000v200_v0%18%266")],
	"Dark flash all magnitude": [("_day6dpfdfall_","_v0D1000v200_v0%18%266"),("_day6dpfdf1_","_v0D1000v200_v0%18%266"),("_day6dpfdf2_","_v0D1000v200_v0%18%266"),("_day6dpfdf3_","_v0D1000v200_v0%18%266")],
	"Dark flash all latency": [("_day6dpfdfall_","_v0D1000v200_v0%18%266"),("_day6dpfdf1_","_v0D1000v200_v0%18%266"),("_day6dpfdf2_","_v0D1000v200_v0%18%266"),("_day6dpfdf3_","_v0D1000v200_v0%18%266")]
	}

stimcombosfull = {
	"Light flash frequency [day & night]": [("_day5dpflfday_", "_v1000D1000v200_"),("_day6dpflfnight_","_v1000D1000v200_"),("_day5dpflfday_", "_v1000_"),("_day6dpflfnight_","_v1000_"),("_day6dpfmslf_", "_v1000D1000v200_")],
	"Light flash magnitude [day & night]": [("_day5dpflfday_", "_v1000D1000v200_"),("_day6dpflfnight_","_v1000D1000v200_"),("_day5dpflfday_", "_v1000_"),("_day6dpflfnight_","_v1000_"),("_day6dpfmslf_", "_v1000D1000v200_")],
	"Light flash latency [day & night]": [("_day5dpflfday_", "_v1000D1000v200_"),("_day6dpflfnight_","_v1000D1000v200_"),("_day5dpflfday_", "_v1000_"),("_day6dpflfnight_","_v1000_"),("_day6dpfmslf_", "_v1000D1000v200_")],
	"Tap habituation frequency [day]": [("_day5dpfhab1_","_a1f1000d5p_"), ("_day5dpfhab2wdf_","_a1f1000d5p_"), ("_day5dpfhab3wdf_","_a1f1000d5p_")],
	"Tap habituation magnitude [day]": [("_day5dpfhab1_","_a1f1000d5p_"), ("_day5dpfhab2wdf_","_a1f1000d5p_"), ("_day5dpfhab3wdf_","_a1f1000d5p_")],
	"Tap habituation latency [day]": [("_day5dpfhab1_","_a1f1000d5p_"), ("_day5dpfhab2wdf_","_a1f1000d5p_"), ("_day5dpfhab3wdf_","_a1f1000d5p_")],
	"Tap post habituation frequency [day]": [("_day5dpfhab3post_","_a1f1000d5p_"), ("_day5dpfhab2post_","_a1f1000d5p_"), ("_day5dpfhab3post_","_a1f1000d5p_")],
	"Tap post habituation magnitude [day]": [("_day5dpfhab3post_","_a1f1000d5p_"), ("_day5dpfhab2post_","_a1f1000d5p_"), ("_day5dpfhab3post_","_a1f1000d5p_")],
	"Tap post habituation latency [day]": [("_day5dpfhab3post_","_a1f1000d5p_"), ("_day5dpfhab2post_","_a1f1000d5p_"), ("_day5dpfhab3post_","_a1f1000d5p_")],
	"Prepulse tap response frequency [day & night]": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a1%89%97"),("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a1%89%97")],
	"Prepulse tap response magnitude [day & night]": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a1%89%97"),("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a1%89%97")],
	"Prepulse tap response latency [day & night]": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a1%89%97"),("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a1%89%97")],
	"Weak tap response frequency [day]": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day5dpfppi_","_a001f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5p_a001%4%12")],
	"Weak tap response magnitude [day]": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day5dpfppi_","_a001f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5p_a001%4%12")],
	"Weak tap response latency [day]": [("_day5dpfppi_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day5dpfppi_","_a001f1000d5p_a001%4%12"),("_day5dpfppi_","_a001f1400d5p_a001%4%12")],
	"Weak tap response frequency [night]": [("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day6dpfppinight_","_a001f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5p_a001%4%12")],
	"Weak tap response magnitude [night]": [("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day6dpfppinight_","_a001f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5p_a001%4%12")],
	"Weak tap response latency [night]": [("_day6dpfppinight_","_a001f1000d5pD300a1f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5pD300a1f1400d5p_a001%4%12"),("_day6dpfppinight_","_a001f1000d5p_a001%4%12"),("_day6dpfppinight_","_a001f1400d5p_a001%4%12")],
	"Strong tap response frequency [day]": [("_day5dpfppi_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a0f1400d5pD300a1f1400d5p_a1%89%97"),("_day5dpfhab1pre_","_a1f1000d5p_a1%4%12")],
	"Strong tap response magnitude [day]": [("_day5dpfppi_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a0f1400d5pD300a1f1400d5p_a1%89%97"),("_day5dpfhab1pre_","_a1f1000d5p_a1%4%12")],
	"Strong tap response latency [day]": [("_day5dpfppi_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day5dpfppi_","_a0f1400d5pD300a1f1400d5p_a1%89%97"),("_day5dpfhab1pre_","_a1f1000d5p_a1%4%12")],
	"Strong tap response frequency [night]": [("_day6dpfppinight_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a0f1400d5pD300a1f1400d5p_a1%89%97")],
	"Strong tap response magnitude [night]": [("_day6dpfppinight_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a0f1400d5pD300a1f1400d5p_a1%89%97")],
	"Strong tap response latency [night]": [("_day6dpfppinight_","_a0f1000d5pD300a1f1000d5p_a1%89%97"),("_day6dpfppinight_","_a0f1400d5pD300a1f1400d5p_a1%89%97")],
	"MS weak tap response L magnitude [day]": [("_day6dpfmslf_","_v1000a001f1000d5pD995v200_a001%4%12")],
	"MS weak tap response NL magnitude [day]": [("_day6dpfmslf_","_v200a001f1000d5pD995v200_a001%4%12")],
	"Dark flash block 1 frequency": [("_day6dpfdf1_","_v0D1000v200_v0%18%266")],
	"Dark flash block 1 magnitude": [("_day6dpfdf1_","_v0D1000v200_v0%18%266")],
	"Dark flash block 1 latency": [("_day6dpfdf1_","_v0D1000v200_v0%18%266")],
	"Dark flash block 2 frequency": [("_day6dpfdf2_","_v0D1000v200_v0%18%266")],
	"Dark flash block 2 magnitude": [("_day6dpfdf2_","_v0D1000v200_v0%18%266")],
	"Dark flash block 2 latency": [("_day6dpfdf2_","_v0D1000v200_v0%18%266")],
	"Dark flash block 3 frequency": [("_day6dpfdf3_","_v0D1000v200_v0%18%266")],
	"Dark flash block 3 magnitude": [("_day6dpfdf3_","_v0D1000v200_v0%18%266")],
	"Dark flash block 3 latency": [("_day6dpfdf3_","_v0D1000v200_v0%18%266")],
	"Dark flash block 1 start frequency": [("_day6dpfdf1a_","_v0D1000v200_v0%18%266")],
	"Dark flash block 1 start magnitude": [("_day6dpfdf1a_","_v0D1000v200_v0%18%266")],
	"Dark flash block 1 start latency": [("_day6dpfdf1a_","_v0D1000v200_v0%18%266")],
	"Dark flash block 1 end frequency": [("_day6dpfdf1b_","_v0D1000v200_v0%18%266")],
	"Dark flash block 1 end magnitude": [("_day6dpfdf1b_","_v0D1000v200_v0%18%266")],
	"Dark flash block 1 end latency": [("_day6dpfdf1b_","_v0D1000v200_v0%18%266")],
	"Dark flash block 2 start frequency": [("_day6dpfdf2a_","_v0D1000v200_v0%18%266")],
	"Dark flash block 2 start magnitude": [("_day6dpfdf2a_","_v0D1000v200_v0%18%266")],
	"Dark flash block 2 start latency": [("_day6dpfdf2a_","_v0D1000v200_v0%18%266")],
	"Dark flash block 2 end frequency": [("_day6dpfdf2b_","_v0D1000v200_v0%18%266")],
	"Dark flash block 2 end magnitude": [("_day6dpfdf2b_","_v0D1000v200_v0%18%266")],
	"Dark flash block 2 end latency": [("_day6dpfdf2b_","_v0D1000v200_v0%18%266")],
	"Dark flash block 3 start frequency": [("_day6dpfdf3a_","_v0D1000v200_v0%18%266")],
	"Dark flash block 3 start magnitude": [("_day6dpfdf3a_","_v0D1000v200_v0%18%266")],
	"Dark flash block 3 start latency": [("_day6dpfdf3a_","_v0D1000v200_v0%18%266")],
	"Dark flash block 3 end frequency": [("_day6dpfdf3b_","_v0D1000v200_v0%18%266")],
	"Dark flash block 3 end magnitude": [("_day6dpfdf3b_","_v0D1000v200_v0%18%266")],
	"Dark flash block 3 end latency": [("_day6dpfdf3b_","_v0D1000v200_v0%18%266")],
	"Dark flash all frequency": [("_day6dpfdfall_","_v0D1000v200_v0%18%266")],
	"Dark flash all magnitude": [("_day6dpfdfall_","_v0D1000v200_v0%18%266")],
	"Dark flash all latency": [("_day6dpfdfall_","_v0D1000v200_v0%18%266")]
	}

direction_swaps = ["rhofrac", "interboutinterval", "fullboutdatamaxloc", "responselatency", "interbouttime", "SLEEP"]
# latency is only going to be swapped when it is being combined with other information, otherwise it is not swapped in the stimulus data

for sfilename in glob.glob("*out"):
#for sfilename in glob.glob("deaf1_c207y_box12_09_22_20_VM_linearmodel_deaf1C207YBox12_het_vs_deaf1C207YBox12_hom_7517615_4294967294.out"):
#for sfilename in glob.glob("*kiaa0232Box12_wt_vs_kiaa0232Box12_hom_11693279_4294967294*out"):
	pvallist = []
	sfile = open(sfilename, 'r')
	#if sfile == "":
	#	sfile = sfilename
	#outfile = open("ssmd" + sfilename + ".csv", 'w')
	lines = sfile.readlines()
	for line in lines:
		#anova:  ribgraph_mean_time_day2dfall_numberofbouts_3600_controlgroup-het.data : N of control, test, Mean of array control, test, Variance of array control, test, SSMD, H-stat, P-value:  43 35 1196.7162790697676 1193.72 272040.68322336394 316867.31017142854 0.0039044380458149457 0.10660666975064714 0.7440409692821573
		if line.startswith("anova:"): # automatically skipping all the failed ones
			pval = line.strip().split()[-1] # I always leave pval at the end
			ssmd = line.strip().split()[-3] # I always leave pval at the end
			graph = line.strip().split(":")[1].strip()
			fulldata = "skipping now"#line.strip().split(":")[3].strip()
			#if float(pval) < cutoff:
			pvallist.append([float(pval), graph, fulldata, ssmd])
		else:
			if line.startswith("ribgraph"):
				ribgraph = line.strip() # saving and then it will match up with the next instance of hitting the "mutornot"
			if line.startswith("linear model failed"):
				ribgraph = line.split(":")[1].strip() # saving and then it will match up with the next instance of hitting the "mutornot"
				pvallist.append([float(1000.0), ribgraph, ""])
			if line.startswith("mutornot[T.wt] "):
				if len(line.split()) > 3:
					lmmpval = line.split()[4]
					coef = line.split()[1]
					if float(lmmpval) == 0:
						lmmpval = 0.001
					#if float(lmmpval) < cutoff:
					pvallist.append([float(lmmpval), ribgraph, ""])
	for stat in pvallist:
		if stat[2] == "": # if it's lmm data
			for s in range(0, len(pvallist)):
				if stat[1] == pvallist[s][1]: # found the lmm data
					pvallist[s][0] = stat[0] # replacing pval in the non-lmm data with lmm pval
	filteredpvallist0 = [x for x in pvallist if not x[2]==""] # eliminate all the original lmms
	#filteredpvallist = [x for x in filteredpvallist0 if not x[0]>cutoff] # filter by the cutoff
	filteredpvallistbase = []
	for final in filteredpvallist0:
		# getting rid of the data that is really short time windows with large binning, as it just adds noise
		if "600" in str(final[1]) and "transition" in str(final[1]):
			continue
		ssmd = -1*float(final[3])
		sig = ""
		if final[0]<0.05:
			sig = "Significant"
		else:
			sig = "Not Significant"
		assaytype = ""
		for x in direction_swaps:
			if x in str(final[1]):
				ssmd = -1*float(ssmd)
		for cats in nonstimcombos.keys():
			graphs = nonstimcombos[cats]
			for g in graphs:
				if "day" in cats:
					if g in str(final[1]) and "day" in str(final[1]):
						assaytype=cats
						continue
			for g in graphs:
				if "night" in cats:
					if g in str(final[1]) and "night" in str(final[1]):
						assaytype=cats
						continue
			if assaytype != "":
				continue
		if assaytype == "":
			continue
		innerlist = [str(final[1]),str(ssmd), assaytype, sig]
		#print(final)
		#print(innerlist)
		filteredpvallistbase.append(innerlist)
#	print(filteredpvallist)
	df = pd.DataFrame(filteredpvallistbase,columns=["Assay", "SSMD","Assay Type","Significance"])
	df.to_csv("ssmd_baseline_" + sfilename + ".csv")

	filteredpvalliststim = []
	for final in filteredpvallist0:
		# getting rid of the data that is really short time windows with large binning, as it just adds noise
		skip = False
		for s in skiplist:
			if s in str(final[1]):
				skip = True
		if skip == True:
			continue
		ssmd = -1*float(final[3])
		sig = ""
		if final[0]<0.05:
			sig = "Significant"
		else:
			sig = "Not Significant"
		assaytype = ""
		for x in direction_swaps:
			if x in str(final[1]):
				ssmd = -1*float(ssmd)
		for cats in stimcombos.keys():
			#print(cats)
			graphs = stimcombos[cats]
			for g in graphs:
			#	print(g)
				if "frequency" in cats:
					#print("TEST", cats)
					#print("TEST2", g[0], g[1], str(final[1]))
					if g[0] in str(final[1]) and g[1] in str(final[1]) and "responsefrequency" in str(final[1]):
						assaytype=cats
						continue
				if "latency" in cats:
					if g[0] in str(final[1]) and g[1] in str(final[1]) and "latency" in str(final[1]):
						assaytype=cats
						continue
				if "magnitude" in cats:
					if g[0] in str(final[1]) and g[1] in str(final[1]) and "latency" not in str(final[1]) and "responsefrequency" not in str(final[1]):
						assaytype=cats
						continue
				else:
					if g[0] in str(final[1]) and g[1] in str(final[1]):
						assaytype=cats
						continue
			if assaytype != "":
				continue
		if assaytype == "":
			continue
		innerlist = [str(final[1]),str(ssmd), assaytype, sig]
#		print(final)
#		print(innerlist)
		filteredpvalliststim.append(innerlist)
#	print(filteredpvalliststim)
	df = pd.DataFrame(filteredpvalliststim,columns=["Assay","SSMD","Assay Type","Significance"])
	df.to_csv("ssmd_stimuli_" + sfilename + ".csv")

	filteredpvalliststimfull = []
	for final in filteredpvallist0:
		# getting rid of the data that is really short time windows with large binning, as it just adds noise
		skip = False
		for s in skiplist:
			if s in str(final[1]):
				skip = True
		if skip == True:
			continue
		ssmd = -1*float(final[3])
		sig = ""
		if final[0]<0.05:
			sig = "Significant"
		else:
			sig = "Not Significant"
		assaytype = ""
	#	for x in direction_swaps:
	#		if x in str(final[1]):
	#			ssmd = -1*float(ssmd)
		for cats in stimcombosfull.keys():
			#print(cats)
			graphs = stimcombosfull[cats]
			for g in graphs:
			#	print(g)
				if "frequency" in cats:
					#print("TEST", cats)
					#print("TEST2", g[0], g[1], str(final[1]))
					if g[0] in str(final[1]) and g[1] in str(final[1]) and "responsefrequency" in str(final[1]):
						assaytype=cats
						continue
				if "latency" in cats:
					if g[0] in str(final[1]) and g[1] in str(final[1]) and "latency" in str(final[1]):
						assaytype=cats
						continue
				if "magnitude" in cats:
					if g[0] in str(final[1]) and g[1] in str(final[1]) and "latency" not in str(final[1]) and "responsefrequency" not in str(final[1]):
						assaytype=cats
						continue
			if assaytype != "":
				continue
		if assaytype == "":
			continue
		innerlist = [str(final[1]),str(ssmd), assaytype, sig]
#		print(final)
#		print(innerlist)
		filteredpvalliststimfull.append(innerlist)
#	print(filteredpvalliststim)
	df = pd.DataFrame(filteredpvalliststimfull,columns=["Assay","SSMD","Assay Type","Significance"])
	df.to_csv("ssmd_stimuli_expanded_" + sfilename + ".csv")


	filteredpvalliststimnorm = []
	for final in filteredpvallist0:
		# getting rid of the data that is really short time windows with large binning, as it just adds noise
		skip = False
		for s in skiplist:
			if s in str(final[1]):
				skip = True
		if skip == True:
			continue
		ssmd = -1*float(final[3])
		sig = ""
		if final[0]<0.05:
			sig = "Significant"
		else:
			sig = "Not Significant"
		assaytype = ""
	#	for x in direction_swaps:
	#		if x in str(final[1]):
	#			ssmd = -1*float(ssmd)
		for cats in stimcombosnorm.keys():
			#print(cats)
			graphs = stimcombosnorm[cats]
			for g in graphs:
			#	print(g)
				if "frequency" in cats:
					#print("TEST", cats)
					#print("TEST2", g[0], g[1], str(final[1]))
					if g[0] in str(final[1]) and g[1] in str(final[1]) and "responsefrequency" in str(final[1]):
						assaytype=cats
						continue
				if "latency" in cats:
					if g[0] in str(final[1]) and g[1] in str(final[1]) and "latency" in str(final[1]):
						assaytype=cats
						continue
				if "magnitude" in cats:
					if g[0] in str(final[1]) and g[1] in str(final[1]) and "latency" not in str(final[1]) and "responsefrequency" not in str(final[1]):
						assaytype=cats
						continue
			if assaytype != "":
				continue
		if assaytype == "":
			continue
		innerlist = [str(final[1]),str(ssmd), assaytype, sig]
#		print(final)
#		print(innerlist)
		filteredpvalliststimnorm.append(innerlist)
#	print(filteredpvalliststim)
	df = pd.DataFrame(filteredpvalliststimnorm,columns=["Assay","SSMD","Assay Type","Significance"])
	df.to_csv("ssmd_stimuli_standard_" + sfilename + ".csv")
