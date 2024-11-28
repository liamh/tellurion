from . import orbit
from . import force
from . import propagate
import astropy.time
import astropy.units as u
import numpy as np

################################################################################
## Cartesian
################################################################################

# Cartesian propagation
ex1 = orbit.new_cart([5740.13268349499, 3314.06715, 0.0,
                -2.75082683526322, 4.7645718414998, 5.50165367052644],
               astropy.time.Time('2022-06-01T12:00:00.000000'),
               force.setgravity(0,0))


# All these are the same as ex1.pvt: PVT(ex1.cart), PVT(ex1.pvt.pvtq), PVT(ex1.pvt.ork), PVT(*ex1.pvt.makenp())
# In [3]: ex1.pvt
# Out[3]: <PVT position: [5740.13268349499, 3314.06715, 0.0] (km) velocity:[-2.75082683526322, 4.7645718414998, 5.50165367052644] (km/s) epoch 2022-06-01T12:00:00.000 (UTC)>
# In [4]: ex1.pvt.pvtq
# Out[4]: <TQuantity ([5740.13268349, 3314.06715   ,    0.        ], [-2.75082684,  4.76457184,  5.50165367]) (km, km / s), time=2022-06-01T12:00:00.000>
# In [8]: ex1.pvt.ork
# Out[8]: <TimeStampedPVCoordinates: {2022-06-01T12:00:00.000, P(5740132.68349499, 3314067.15, 0.0), V(-2750.82683526322, 4764.5718414998, 5501.65367052644), A(0.0, 0.0, 0.0)}>
# elementval(ex1.cart, "sma")
# Out[20]: <Quantity 6672.37441023 km>

# ex1.cartesian()

def propdemo():
    """Propagation demo on ex1
    After running this demo, ex1.prop.orbit, ex1.prop.pvt, and
    ex1.prop.ephem will have variants of the the 10-step ephemeris.

    import astrodynamics.example as exmp
    exmp.propdemo()
    exmp.ex1.prop.ephem
    <TimeSeries length=10>
            time               position               velocity
                                  km                   km / s
            Time              float64[3]             float64[3]
    ------------------- ---------------------- ----------------------
    2022-06-01 12:02:00  5354.637 ..   658.032 -3.663509 ..  5.447524
    2022-06-01 12:02:15  5298.863 ..   739.639 -3.772863 ..  5.433178
    2022-06-01 12:02:30  5241.457 ..   821.019 -3.881049 ..  5.417159
    2022-06-01 12:02:45  5182.438 ..   902.146 -3.988033 ..  5.399474
    2022-06-01 12:03:00  5121.822 ..   982.995 -4.093782 ..  5.380127
    2022-06-01 12:05:00  4581.815 ..  1616.826 -4.891449 ..  5.166423
    2022-06-01 12:10:00  2865.500 ..  3036.847 -6.432387 ..  4.203822
    2022-06-01 12:20:00 -1359.960 ..  4646.558 -7.074062 ..  0.948420
    2022-06-01 12:25:00 -3358.332 ..  4647.312 -6.115174 .. -0.941270
    2022-06-01 12:30:00 -4956.792 ..  4094.286 -4.436284 .. -2.706752

    From ephem, use makenp() to turn into numpy arrays,
    ex1.prop.ephem[5:10].makenp()
    propagate.pvts(ex1.prop.ephem) # same as ex1.prop.pvt

    This function will show a warning which should be ignored.
    """
    propagate.prop(ex1,300.0,maximum_tof=86400.0)
    propagate.prop(ex1,[1200.0,1500.0,1800.0])
    propagate.prop(ex1,600.0)
    propagate.prop(ex1, np.linspace(2.0, 3.0, num=5, endpoint=True)*u.minute)

################################################################################
## Kepler elements
################################################################################

# Build and convert a Kepler
ex2 = orbit.new_kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "raan":217.4, "ma":7.25},
                 astropy.time.Time('2023-09-14T08:30:00'),
                 force.setgravity(0, 0))

## Example
# prop(o.ex2, np.linspace(0, 24, num=5)*u.hour) # This gives a warning, can be ignored
# o.ex2.prop.ephem
# <TimeSeries length=5>
#         time               position               velocity
#                               km                   km / s
#         Time              float64[3]             float64[3]
# ------------------- ---------------------- ----------------------
# 2023-09-14 08:30:00 -4100.052 .. -4764.959 -5.636147 ..  0.734345
# 2023-09-14 14:30:00 -5300.257 .. -4452.057 -4.456035 ..  1.892044
# 2023-09-14 20:30:00 -6192.600 .. -3879.939 -3.052240 ..  2.910624
# 2023-09-15 02:30:00 -6737.221 .. -3089.289 -1.537353 ..  3.728167
# 2023-09-15 08:30:00 -6920.676 .. -2132.309 -0.019043 ..  4.309723
#
# tselements(o.ex2.prop, ["sma","ecc","inc"])
# <TimeSeries length=5>
#         time               sma                ecc                inc
#                             km                                   deg
#         Time             float64            float64            float64
# ------------------- ----------------- ------------------- ------------------
# 2023-09-14 08:30:00 8000.000000000001 0.10000000000000005               42.0
# 2023-09-14 14:30:00 7999.984739653951 0.09999739630636867  41.99999999999989
# 2023-09-14 20:30:00 7999.975777150673 0.09999489374642954 41.999999999999865
# 2023-09-15 02:30:00 7999.961193141438 0.09999236890210698  42.00000000000003
# 2023-09-15 08:30:00 7999.949375719636  0.0999895810073231 42.000000000000306
#
# tselements(o.ex2.prop, ["altper","altapo"])
# <TimeSeries length=5>
#         time              altper            altapo
#                             km                km
#         Time             float64           float64
# ------------------- ----------------- ------------------
# 2023-09-14 08:30:00 821.8635400000001  2421.863540000002
# 2023-09-14 14:30:00 821.8706351978723  2421.825924110028
# 2023-09-14 20:30:00  821.882589340481  2421.796044960863
# 2023-09-15 02:30:00 821.8896623142995  2421.759803968576
# 2023-09-15 08:30:00  821.901329561633 2421.7245018776366
