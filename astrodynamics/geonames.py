"""
Find lat, lon, elev of places by name; see
https://gitlab.com/-/snippets/3743091 for info on how to get a
geonames id.
 from astrodynamics import geonames
 utspaceeng = geonames.location('speedway parking garage')
"""

import geocoder
import requests

userid = ''
def gnuserid():
    global userid
    if userid=='':
        userid = input("Enter your geonames userid: ")
    return userid

def location(place):
    resp = geocoder.geonames(place, key=gnuserid())
    props = resp.geojson['features'][0]['properties']
    lat = float(props['lat'])
    lon = float(props['lng'])
    return [lat, lon, elev(lat, lon), props['address'], props['state'], props['country']]

def elev(lat, lon):
    r = requests.get(f"https://api.opentopodata.org/v1/aster30m?locations={lat},{lon}")
    data = r.json()
    return data['results'][0]['elevation']
