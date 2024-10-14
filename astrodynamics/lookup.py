"""
Find lat, lon, elev of places by name; see
https://gitlab.com/-/snippets/3743091 for info on how to get a
geonames id.
 import lookup
 lookup.geonames_userid = 'myusername'
 utspaceeng = lookup.location('speedway parking garage')
"""

import geocoder
import requests
def location(place):
    global geonames_userid
    resp = geocoder.geonames(place, key=geonames_userid)
    props = resp.geojson['features'][0]['properties']
    lat = float(props['lat'])
    lon = float(props['lng'])
    return [lat, lon, elev(lat, lon), props['address'], props['state'], props['country']]

def elev(lat, lon):
    r = requests.get(f"https://api.opentopodata.org/v1/aster30m?locations={lat},{lon}")
    data = r.json()
    return data['results'][0]['elevation']
