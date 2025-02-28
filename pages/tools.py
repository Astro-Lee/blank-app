import streamlit as st

#!/usr/bin/env python
from datetime import datetime as dt
from datetime import timedelta as td
import astropy.units as u
from astropy.coordinates import SkyCoord


def time_diff(trigger: str, time: str):
    trigger = dt.fromisoformat(trigger)
    delta = (dt.fromisoformat(time) - trigger).total_seconds()

    st.write(f'{delta:.0f} seconds')
    st.write(f'{delta/60:.2f} minutes')
    st.write(f'{delta/3600:.2f} hours')
    st.write(f'{delta/3600/24:.2f} days')

def separation(coord1: str, coord2: str, deg=False):
    if deg:
        coord1 = SkyCoord(coord1, unit=(u.deg,u.deg))
        coord2 = SkyCoord(coord2, unit=(u.deg,u.deg))
    else:
        coord1 = SkyCoord(coord1, unit=(u.hourangle, u.deg))
        coord2 = SkyCoord(coord2, unit=(u.hourangle, u.deg))
    
    sep = coord1.separation(coord2)
    
    st.write(f'{sep.to(u.arcsec):.2f}')
    st.write(f'{sep.to(u.arcmin):.2f}')
    st.write(f'{sep:.2f}')
    
def flux_error(mean_index=-10.8706, lower_index=-10.9855,upper_index=-10.7732):
    mean_flux = 10 ** mean_index
    upper_flux = 10 ** upper_index
    lower_flux = 10 ** lower_index

    st.write(f'{mean_flux:.2e} (+{upper_flux - mean_flux:.2e}, {lower_flux-mean_flux:.2e}) erg/cm^2/s')

st.title('Tools')

trigger_time = st.text_input("Trigger Time", value="2024-08-25T00:00:00")
time = st.text_input("Time", value="2024-08-25T00:00:00")
st.button('Calculate', key='time_diff',on_click=time_diff(trigger_time, time))

coord1 = st.text_input("Coordinate 1", value="00h00m00s +00d00m00s")
coord2 = st.text_input("Coordinate 2", value="00h00m00s +00d00m00s")
in_deg = st.checkbox('Coordinates in degrees', key='deg')
st.button('Calculate', key='separation',on_click=separation(coord1,coord2,in_deg))


mean_index = st.number_input("Mean Index", value=-10.8706)
lower_index = st.number_input("Lower Index", value=-10.9855)
upper_index = st.number_input("Upper Index", value=-10.7732)

st.button('Calculate', key='flux_error', on_click=flux_error(mean_index, lower_index,upper_index))