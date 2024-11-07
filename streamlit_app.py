import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import smplotlib

import requests
import json
# from IPython.display import IFrame, Image
from astropy.coordinates import SkyCoord
import astropy.units as u
import urllib.parse
from datetime import datetime, timezone

st.title("GMG Follow-Up Plan")

# Sidebar inputs for source information
st.sidebar.text_input("Source Name", key="source_name", value="GRB 240825A")
st.sidebar.text_input("T0", key="T0", value="2024-08-25T15:52:59")

# RA, DEC, error radius, and unit in a single row with four columns
col1, col2, col3, col4 = st.sidebar.columns(4)
col1.text_input("RA", key="RA", value="344.57192")
col2.text_input("DEC", key="DEC", value="1.02686")
col3.text_input("Err. Radius", key="Err", value="0.74")
col4.radio("Unit", ["arcmin", "arcsec", "deg"], key="unit")

# Retrieve session state values
source_name = st.session_state.source_name
T0 = st.session_state.T0
RA = st.session_state.RA
DEC = st.session_state.DEC
Err = st.session_state.Err
unit = st.session_state.unit

# Display the input values in the main page
# st.write("Source Name:", source_name)
# st.write("T0:", T0)
# st.write("RA:", RA)
# st.write("DEC:", DEC)
# st.write(f"Err. Radius: {Err} {unit}")

coords = SkyCoord(RA, DEC, unit=(u.deg, u.deg), frame='icrs', obstime=T0)

Err_deg = (Err*u.arcmin).to(u.deg)

raDeg = coords.ra.deg
decDeg = coords.dec.deg

lDeg = coords.galactic.l.deg
bDeg = coords.galactic.b.deg

raJ2000 = coords.ra.to_string(u.hourangle,pad=True,sep=' ',precision=2)
decJ2000 = coords.dec.to_string(u.deg,pad=True,alwayssign=True,sep=' ',precision=1)
coordsJ2000 = f'{raJ2000} {decJ2000}'

# st.write(f"""{source_name}
         
# T0 at {T0}

# Localization (J2000):

# RA  =  {raJ2000}

# DEC = {decJ2000}

# Err = {Err} {unit}
# """)

from bs4 import BeautifulSoup
params = {
    'Coords': coordsJ2000,
    'equinox': '2000',
    'ascii': '1' # 1 for ascii, 0 for html
}
response = requests.post('https://www.swift.ac.uk/analysis/nhtot/donhtot.php', params=params)

# Check if the request was successful
if response.status_code == 200:
    # Save the content as a GIF file
    # with open('output.gif', 'wb') as file:
    #     file.write(response.content)
    # print("Image downloaded and saved as 'output.gif'.")
    nH_url = response.url
else:
    print(f"Failed to retrieve the image. Status code: {response.status_code}")
    
# print(nH_url)

# 使用 BeautifulSoup 解析 HTML
soup = BeautifulSoup(response.text, 'html.parser')

# 查找 <pre> 标签中的内容
pre_tag = soup.find('pre')

# 获取 <pre> 标签中的文本内容
if pre_tag:
    absorption = pre_tag.get_text()
    # st.write(absorption)
else:
    st.write("No <pre> tag found.")

# Define the URL
url = """https://aladin.cds.unistra.fr/AladinLite/?fov=0.10&survey=CDS/P/PanSTARRS/DR1/color-z-zg-g&overlays=[{ "type": "graphic", "name": "uncertainty", "color": "red", "lineWidth": 4, "items": [ { "type": "circle", "ra": %s, "dec": %s, "radius": %s } ]}]&target=%s"""%(raDeg,decDeg,Err_deg,coordsJ2000)

encoded_url = urllib.parse.quote(url, safe=":/?=&")
# st.write(encoded_url)

# Get the current UTC time
current_utc_time = datetime.now(timezone.utc)
current_year = current_utc_time.year
current_month = current_utc_time.month
current_day = current_utc_time.day

# Query parameters as a dictionary
params = {
    'action': 'showImage',
    'form[mode]': '1',
    'form[day]': f'{current_day}',
    'form[month]': f'{current_month}',
    'form[year]': f'{current_year}',
    'form[obs_name]': 'Lijiang Observatory (China)',
    'form[coordlist]': f'{coordsJ2000}',
    'form[coordfile]': '',
    'form[paramdist]': '2',
    'form[format]': 'gif',
    'submit': 'Retrieve'
}

# Send a GET request
response = requests.get('https://astro.ing.iac.es/staralt/index.php', params=params)

# Check if the request was successful
if response.status_code == 200:
    # Save the content as a GIF file
    # with open('output.gif', 'wb') as file:
    #     file.write(response.content)
    # print("Image downloaded and saved as 'output.gif'.")
    staralt_url = response.url
else:
    print(f"Failed to retrieve the image. Status code: {response.status_code}")
    
# st.write(staralt_url)

st.write(
        f'<iframe width="850px" height="680px" src="{staralt_url}"></iframe>'
        '<iframe width="1400" height="900" src="http://weather.gmg.org.cn:9000">',
        unsafe_allow_html=True,
    )

st.write(f"""{source_name}

T0 at {T0}

Localization (J2000):

RA  =  {raJ2000}

DEC = {decJ2000}

Err = {Err} {unit}

Galactic l, b = ({lDeg:.1f}, {bDeg:.1f}) deg

{absorption}

Visibility: {staralt_url}

Finding Chart: {encoded_url}

Weather: http://weather.gmg.org.cn:9000
""")

import requests as req
import tempfile
import xml.etree.ElementTree as ET
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
from astropy.io import fits
from datetime import datetime, timedelta
from tabulate import tabulate


target_name = source_name
circular_num = 37274
TA_last_name = 'Gupta'
T0 = T0

results = ['https://stdweb.astro0871.top/tasks/2/download/candidates/J225817.27+010136.8.cutout',
          ]

obs_time = []
req_time = []
ra = []
dec = []
filter = []
mag = []
mag_err = []
mag_limit = []

for result in results:
    rsp = req.get(result)
    if rsp.status_code == 200:
        
        with tempfile.NamedTemporaryFile() as tmpf:
            tmpf.write(rsp.content)
            
            with fits.open(tmpf,mode='update') as hdu:
                #hdu.info()
                target = hdu[0].header
                obs_info = hdu[1].header
                
    coord = SkyCoord(target['RA'], target['DEC'], unit=(u.deg, u.deg), obstime=obs_info['DATE-OBS'])
    RA = coord.ra.to_string(u.hourangle,pad=True,sep=' ',precision=2)
    DEC = coord.dec.to_string(u.deg,pad=True,sep=' ',precision=1,alwayssign=True)
    
    obs_time.append(obs_info['DATE-OBS'])
    req_time.append(obs_info['REQTIME'])
    ra.append(RA)
    dec.append(DEC)
    filter.append(target['mag_filter_name '][0])
    mag.append(target['mag_calib'])
    mag_err.append(target['mag_calib_err'])
    mag_limit.append(target['mag_limit'])

info = pd.DataFrame(data={
    'obs_time': obs_time,
    'req_time': req_time,
    'ra': ra,
    'dec': dec,
    'filter': filter,
    'mag': mag,
    'mag_err': mag_err,
    'mag_limit': mag_limit}).sort_values(by='obs_time', ignore_index=True)

info['mag'] = info['mag'].round(2)

info['mag_err'] = info['mag_err'].round(2).apply(lambda x: 0.01 if x < 0.01 else x)

info['mag_limit'] = info['mag_limit'].round(1)

start_obs = f"{(datetime.fromisoformat(info['obs_time'][0]) - datetime.fromisoformat(T0)).total_seconds() / 3600:.2f}"

for i, row in info.iterrows():
    
    info.loc[i, 'hour_diff'] = f"{(datetime.fromisoformat(row['obs_time']) + timedelta(seconds=float(row['req_time'])/2) - datetime.fromisoformat(T0)).total_seconds() / 3600:.2f}"

    info.loc[i,'mag_fmt'] = f"{row['mag']} ± {row['mag_err']}"
info

# st.write(info)


rsp = req.get(url='https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust', 
              params={'locstr': coord.to_string(style='hmsdms', precision=2, alwayssign=True)[1:],
                      'regSize': '2.0'})

if rsp.status_code == 200:
#    print('SF11 Reddening map 请求成功！')
    with tempfile.NamedTemporaryFile() as fp:
        fp.write(rsp.content)
        fp.seek(0)
        tree = ET.parse(fp)

reddening = f'{float(tree.find(".//refPixelValueSandF").text.strip().split()[0]):.4f}'

table = tabulate(info[['hour_diff','req_time','filter','mag_fmt','mag_limit']], 
                 headers=['Tmid-T0 [hr]','Exp. [s]','Filter','Mag','5-sigma U.L.'], 
                 tablefmt='grid',numalign='center',stralign='center',showindex=False)

st.write(f"""
{target_name}: GMG Optical Observation
=====================================================================
R.-Z. Li, B.-T. Wang, F.-F. Song, J. Mao,       and J.-M. Bai (YNAO, CAS) report:

We observed the field of {target_name} ({TA_last_name} et al., GCN {circular_num}, T0 at {T0}) using the GMG-2.4m telescope at the Lijiang Observatory. The observation began at {info['obs_time'][0].split('.')[0]}, about {start_obs} hours after the trigger.

The optical counterpart of {target_name}, not visible in the Pan-STARRS1 r-band image, was (marginally/clearly) detected at the coordinates (J2000):
RA  =  {info['ra'][0]}
DEC = {info['dec'][0]}
, with a positional uncertainty of 0.5" or better. The position is consistent with previous results (GCN XXX).

The preliminary analysis results are shown as follows:
{table} 
The given magnitudes are derived based on calibration against Pan-STARRS1 field stars, and are not corrected for the expected Galactic foreground extinction, corresponding to a reddening of E(B-V) = {reddening} mag in the direction of the optical counterpart (Schlafly & Finkbeiner 2011).
""")

table = tabulate(info[['hour_diff','req_time','filter','mag_limit']], 
                 headers=['Tmid-T0 [hr]','Exp. [s]','Filter','5-sigma U.L.'], 
                 tablefmt='grid',numalign='center',stralign='center',showindex=False)

st.write(f"""
{target_name}: GMG Optical Upper Limit
=====================================================================
R.-Z. Li, B.-T. Wang, F.-F. Song, J. Mao,       and J.-M. Bai (YNAO, CAS) report:

We observed the field of {target_name} ({TA_last_name} et al., GCN {circular_num}, T0 at {T0}) using the GMG-2.4m telescope at the Lijiang Observatory. The observation began at {info['obs_time'][0].split('.')[0]}, about {start_obs} hours after the trigger.

No new uncataloged optical source was detected within the XRT/EP-WXT/EP-FXT/VT error circle (GCN XXX).

The preliminary analysis results are shown as follows:
{table}
The given magnitudes are derived based on calibrating against Pan-STARRS1 field stars.
""")

# coords = SkyCoord(RA, DEC, unit=(u.deg, u.deg), frame='icrs', obstime=T0)

# Err_deg = Err.to(u.deg).value

# raDeg = coords.ra.deg
# decDeg = coords.dec.deg

# lDeg = coords.galactic.l.deg
# bDeg = coords.galactic.b.deg

# raJ2000 = coords.ra.to_string(u.hourangle,pad=True,sep=' ',precision=2)
# decJ2000 = coords.dec.to_string(u.deg,pad=True,alwayssign=True,sep=' ',precision=1)
# coordsJ2000 = f'{raJ2000} {decJ2000}'
# print(coordsJ2000)


# add_selectbox = st.sidebar.selectbox(
#     'How would you like to be contacted?',
#     ('arcsec', 'arcmin')
# )








# fig, ax = plt.subplot_mosaic([['a']
#                               ])

# xx = np.linspace(0,2*np.pi,100)
# ax['a'].plot(xx, np.cos(xx))

# st.write(fig)

