"""
Find lat, lon, elev of places by name; see
https://gitlab.com/-/snippets/3743091 for info on how to get a
geonames id.
 from astrodynamics import geonames
 utspaceeng = geonames.location('speedway parking garage')
"""

import geocoder
import requests
import os

userid = ''
def gnuserid():
    '''The Geonames userid, which is necessary to use the Geonames
    placename server and obtained by going to
    https://www.geonames.org/manageaccount and enabling API use with
    box at bottom. Set environment variable in shell export GEONAMES_USERID=xxxx.'''
    global userid
    if userid is None or userid=='':
        userid = os.getenv("GEONAMES_USERID")
    if userid is None or userid=='':
        userid = input("Enter your geonames userid: ")
        print("Avoid prompting with `export GEONAMES_USERID=%s` from shell" % userid)
    return userid

def location(place):
    '''Find the location in Geonames from the name.'''
    resp = geocoder.geonames(place, key=gnuserid())
    props = resp.geojson['features'][0]['properties']
    lat = float(props['lat'])
    lon = float(props['lng'])
    return [lat, lon, elev(lat, lon), props['address'], props['state'], props['country']]

def elev(lat, lon):
    '''Find the elevation from the latitude and longitude from opentopodata.org.'''
    r = requests.get(f"https://api.opentopodata.org/v1/aster30m?locations={lat},{lon}")
    data = r.json()
    return data['results'][0]['elevation']
