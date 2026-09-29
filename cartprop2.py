# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Propagator Comparison Tutorial
#
# This notebook demonstrates the four available propagators: Keplerian, Brouwer-Lyddane, DSST, and Numerical.
# It compares their speed, accuracy, and suitability for different orbital scenarios.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Setup

# %% editable=true slideshow={"slide_type": ""}
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import time
from munch import Munch

import tellurion as tell


# %% [markdown] editable=true slideshow={"slide_type": ""}
# Define test epoch and propagation time arrays for various durations

# %% editable=true slideshow={"slide_type": ""}
newyear = tell.abstime('2025-01-01T00:00:00')

# Short propagation: 1 hour in 5-minute steps
prop_short = np.linspace(5.0*u.minute, 60.0*u.minute, 12)

# Medium propagation: 1 day in 2-hour steps
prop_medium = np.linspace(2.0*u.hour, 24.0*u.hour, 12)

# Long propagation: 10 days in 1-day steps
prop_long = np.linspace(1.0*u.day, 10.0*u.day, 10)


# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Orbital Test Cases
#
# Define three representative orbits: LEO, elliptical, and GEO

# %%
def keppvt(dictels, time):
    """From the dictionary of element values, compute the Kepler
    elements (sma, ecc, inc, argper, raan, and ma or ta), the PVT
    transformation of those elements, and the spherical coordinates."""
    d = {'els': dictels}
    d['kep'] = tell.kepler(tell.allplane(d['els']), time)
    d['pvt'] = d['kep'].pvt()
    d['pvt'].spherical
    return d

leo2 = keppvt({"altper": 525.0*u.km, "altapo": 555.0*u.km,
                       "inc":28.5*u.deg, "argper": 40.0*u.deg,
                       "raan": 40.0*u.deg, "ma": -30.0*u.deg},
                      tell.abstime('2026-01-01 09:30:00'))

vang1 = keppvt({"altper": 600.0*u.km, "altapo": 12000.0*u.km,
                   "inc":36.0*u.deg, "argper": 140.0*u.deg,
                   "raan": 0.0*u.deg, "ma": 100.0*u.deg},
                  tell.abstime('2026-01-01 14:45:00'))

geo1 = keppvt({"memo":1.0*u.rev/u.sday,
                   "ecc":0.0,
                   "inc":0.0*u.deg, "argper": 120.0*u.deg,
                   "raan": 0.0*u.deg, "ma": 0.0*u.deg},
                  tell.abstime('2026-01-01 20:30:00'))

print("LEO Orbit:")
print(f"  sma: {leo2['kep'].elements['sma']:.1f}")
print(f"  eccentricity: {leo2['kep'].elements['ecc']:.4f}")
print(f"  inclination: {leo2['kep'].elements['inc']:.1f}")


# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Propagator 1: Keplerian (Two-Body)
#
# The fastest propagator with no perturbations. Baseline for comparison.

# %% editable=true slideshow={"slide_type": ""}
# Propagate LEO orbit with Keplerian
t0 = time.time()
gen_kep = tell.prepare(leo2['pvt'], 1.0*u.day, propagator='keplerian')
ephem_kep = tell.propagate(gen_kep, prop_medium)
t_kep = time.time() - t0

print(f"Keplerian propagation time: {t_kep*1000:.2f} ms")
print(f"Number of points: {len(ephem_kep)}")
print(f"Position at t=12h: {ephem_kep[6]['position']}")


# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Propagator 2: Brouwer-Lyddane (J2-J5 Perturbations)
#
# Analytic propagator with mean element theory. Much more accurate than Keplerian for perturbed orbits,
# but still very fast. Limited to J2-J5 gravity model. For unknown reasons, this fails on `leo2`, but it works for `vang1`.

# %%
# Propagate LEO orbit with Brouwer-Lyddane (J2-J5)
forceenv_bl = tell.setgravity(5, 0)  # J2-J5 (zonal only)

t0 = time.time()

# Unfortunately, this example results in a non-convergence error
# org.orekit.errors.OrekitException: org.orekit.errors.OrekitException: unable to compute Brouwer-Lyddane mean parameters after 501 iterations
gen_bl = tell.prepare(leo2['pvt'], 1.0*u.day, forceenv=forceenv_bl, propagator='brouwer-lyddane')
ephem_bl = tell.propagate(gen_bl, prop_medium)
t_bl = time.time() - t0

print(f"Brouwer-Lyddane propagation time: {t_bl*1000:.2f} ms")
print(f"Position at t=12h: {ephem_bl[6]['position']}")

# Compare with Keplerian
pos_diff_bl = tell.magdiff(ephem_kep['position'], ephem_bl['position'])
print(f"\nMax position difference vs Keplerian: {np.max(pos_diff_bl):.3f}")
print(f"Mean position difference vs Keplerian: {np.mean(pos_diff_bl):.3f}")


# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Propagator 3: DSST (Semi-Analytical)
#
# Semi-analytical propagator combining numerical integration with perturbation averaging.
# Excellent for long-term propagation with full gravity model. 2-10x faster than numerical integration.

# %%
# Propagate LEO orbit with DSST (full gravity model)
forceenv_dsst = tell.setgravity(20, 20)  # Full 20x20 gravity model

t0 = time.time()
gen_dsst = tell.prepare(leo2['pvt'], 1.0*u.day, forceenv=forceenv_dsst, propagator='dsst')
ephem_dsst = tell.propagate(gen_dsst, prop_medium)
t_dsst = time.time() - t0

print(f"DSST propagation time: {t_dsst*1000:.2f} ms")
print(f"Position at t=12h: {ephem_dsst[6]['position']}")

# Compare with Keplerian
pos_diff_dsst = tell.magdiff(ephem_kep['position'], ephem_dsst['position'])
print(f"\nMax position difference vs Keplerian: {np.max(pos_diff_dsst):.3f}")
print(f"Mean position difference vs Keplerian: {np.mean(pos_diff_dsst):.3f}")


# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Propagator 4: Numerical Integration
#
# Full numerical integration with highest accuracy. Slowest method but best for short-term, high-precision propagation.

# %%
# Propagate LEO orbit with Numerical integration (full gravity model)
forceenv_num = tell.setgravity(20, 20)  # Full 20x20 gravity model

t0 = time.time()
gen_num = tell.prepare(leo2['pvt'], 1.0*u.day, forceenv=forceenv_num, propagator='numerical')
ephem_num = tell.propagate(gen_num, prop_medium)
t_num = time.time() - t0

print(f"Numerical propagation time: {t_num*1000:.2f} ms")
print(f"Position at t=12h: {ephem_num[6]['position']}")

# Compare with Keplerian
pos_diff_num = tell.magdiff(ephem_kep['position'], ephem_num['position'])
print(f"\nMax position difference vs Keplerian: {np.max(pos_diff_num):.3f}")
print(f"Mean position difference vs Keplerian: {np.mean(pos_diff_num):.3f}")


# %% [markdown]
# ## Speed Comparison
#
# Compare execution times and speed ratios

# %%
# Create comparison table
times = {
    'Keplerian': t_kep,
    'Brouwer-Lyddane': t_bl,
    'DSST': t_dsst,
    'Numerical': t_num
}

print("Speed Comparison (LEO, 1 day propagation):")
print("-" * 50)
for name, t in sorted(times.items(), key=lambda x: x[1]):
    ratio = t / t_kep
    print(f"{name:20s}: {t*1000:8.2f} ms  ({ratio:6.1f}x Keplerian)")

print(f"\nDSST speedup over Numerical: {t_num/t_dsst:.1f}x")


# %% [markdown]
# ## Accuracy Comparison
#
# Compare position differences against Numerical (reference) and Keplerian (baseline)

# %%
# Compare all propagators against numerical as reference
print("Position Error vs Numerical Integration (Reference):")
print("-" * 50)

diff_kep_num = tell.magdiff(ephem_kep['position'], ephem_num['position'])
diff_bl_num = tell.magdiff(ephem_bl['position'], ephem_num['position'])
diff_dsst_num = tell.magdiff(ephem_dsst['position'], ephem_num['position'])

errors = {
    'Keplerian': diff_kep_num,
    'Brouwer-Lyddane': diff_bl_num,
    'DSST': diff_dsst_num
}

for name, diffs in errors.items():
    print(f"{name:20s}:")
    print(f"  Max error:  {np.max(diffs):10.3f}")
    print(f"  Mean error: {np.mean(diffs):10.3f}")
    print(f"  Std dev:    {np.std(diffs):10.3f}")
    print()


# %% [markdown]
# ## Long-Term Propagation: DSST Advantage
#
# Demonstrate DSST's advantage over numerical for longer propagations (10 days)

# %%
print("Long-term propagation (10 days):")
print("-" * 50)

# DSST for 10 days
t0 = time.time()
gen_dsst_long = tell.prepare(leo2['pvt'], 10.0*u.day, forceenv=forceenv_dsst, propagator='dsst')
ephem_dsst_long = tell.propagate(gen_dsst_long, prop_long)
t_dsst_long = time.time() - t0

print(f"DSST (10 days):     {t_dsst_long*1000:8.2f} ms")

# Numerical for 10 days
t0 = time.time()
gen_num_long = tell.prepare(leo2['pvt'], 10.0*u.day, forceenv=forceenv_num, propagator='numerical')
ephem_num_long = tell.propagate(gen_num_long, prop_long)
t_num_long = time.time() - t0

print(f"Numerical (10 days): {t_num_long*1000:8.2f} ms")
print(f"\nSpeedup: {t_num_long/t_dsst_long:.1f}x")


# %% [markdown]
# ## Orbital Perturbation Effects
#
# Examine how J2 perturbations affect orbital elements over time

# %%
# Extract RAAN (right ascension of ascending node) at each time step
# This shows the J2 precession effect

raan_kep = []
raan_bl = []
raan_dsst = []
raan_num = []

for i in range(len(ephem_kep)):
    kep_kep = ephem_kep[i].pvt().kepler().elements.to_dict()
    kep_bl = ephem_bl[i].pvt().kepler().elements.to_dict()
    kep_dsst = ephem_dsst[i].pvt().kepler().elements.to_dict()
    kep_num = ephem_num[i].pvt().kepler().elements.to_dict()
    
    raan_kep.append(kep_kep['raan'].value)
    raan_bl.append(kep_bl['raan'].value)
    raan_dsst.append(kep_dsst['raan'].value)
    raan_num.append(kep_num['raan'].value)

times_hours = (prop_medium.to(u.hour).value)

# Plot RAAN evolution
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(times_hours, raan_kep, 'o-', label='Keplerian', linewidth=2)
ax.plot(times_hours, raan_bl, 's-', label='Brouwer-Lyddane', linewidth=2)
ax.plot(times_hours, raan_dsst, '^-', label='DSST', linewidth=2)
ax.plot(times_hours, raan_num, 'd-', label='Numerical', linewidth=2)
ax.set_xlabel('Time (hours)', fontsize=12)
ax.set_ylabel('RAAN (degrees)', fontsize=12)
ax.set_title('RAAN Evolution: Effect of J2 Perturbations (LEO)', fontsize=14)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print("RAAN change over 24 hours:")
print(f"  Keplerian:          {raan_kep[-1] - raan_kep[0]:+.4f}°")
print(f"  Brouwer-Lyddane:    {raan_bl[-1] - raan_bl[0]:+.4f}°")
print(f"  DSST:               {raan_dsst[-1] - raan_dsst[0]:+.4f}°")
print(f"  Numerical:          {raan_num[-1] - raan_num[0]:+.4f}°")


# %% [markdown]
# ## Multi-Orbit Comparison
#
# Propagate different orbit types and compare propagator performance

# %%
# Test all propagators on all orbit types
orbits = {'LEO': leo2, 'Elliptical': vang1, 'GEO': geo1}
propagators = ['keplerian', 'brouwer-lyddane', 'dsst', 'numerical']
forceenvs = {
    'keplerian': None,
    'brouwer-lyddane': tell.setgravity(5, 0),
    'dsst': tell.setgravity(20, 20),
    'numerical': tell.setgravity(20, 20)
}

results = Munch()

for orbit_name, orbit_data in orbits.items():
    results[orbit_name] = Munch()
    
    for prop_type in propagators:
        forceenv = forceenvs[prop_type]
        
        try:
            t0 = time.time()
            gen = tell.prepare(orbit_data['pvt'], 1.0*u.day, 
                              forceenv=forceenv, propagator=prop_type)
            ephem = tell.propagate(gen, prop_medium)
            t_exec = time.time() - t0
            
            results[orbit_name][prop_type] = {
                'time': t_exec,
                'success': True
            }
        except Exception as e:
            results[orbit_name][prop_type] = {
                'time': None,
                'success': False,
                'error': str(e)
            }

# Print results table
print("Propagation Times (ms) across orbit types:")
print("=" * 70)
print(f"{'Orbit Type':<15} {'Keplerian':<15} {'Brouwer-Lyddane':<18} {'DSST':<12} {'Numerical':<12}")
print("-" * 70)

for orbit_name in orbits.keys():
    row = [orbit_name]
    for prop_type in propagators:
        if results[orbit_name][prop_type]['success']:
            t_ms = results[orbit_name][prop_type]['time'] * 1000
            row.append(f"{t_ms:6.2f}")
        else:
            row.append("FAIL")
    print(f"{row[0]:<15} {row[1]:<15} {row[2]:<18} {row[3]:<12} {row[4]:<12}")


# %% [markdown]
# ## Propagator Selection Guide
#
# Summary and recommendations for choosing a propagator

# %%
print("""
PROPAGATOR SELECTION GUIDE
===========================

1. KEPLERIAN (Two-Body)
   Use when:
   - Speed is critical (fastest option)
   - No perturbations needed
   - Quick estimates acceptable
   Speed: ★★★★★ (baseline)
   Accuracy: ★★☆☆☆ (basic)
   Duration: Any

2. BROUWER-LYDDANE (J2-J5)
   Use when:
   - Mean elements (osculating->mean conversion) needed
   - J2-J5 perturbations sufficient (low-order zonal)
   - Medium accuracy required
   - Fast computation essential
   Speed: ★★★★☆ (~same as Keplerian)
   Accuracy: ★★★☆☆ (medium)
   Duration: Days

3. DSST (Semi-Analytical)
   Use when:
   - Long-term propagation (days to weeks)
   - Full gravity model needed (tesseral harmonics)
   - Speed/accuracy tradeoff important
   - 2-10x faster than numerical needed
   Speed: ★★★☆☆ (2-10x slower than Keplerian, 2-10x faster than Numerical)
   Accuracy: ★★★★☆ (high)
   Duration: Weeks-months

4. NUMERICAL INTEGRATION
   Use when:
   - Highest accuracy required
   - Short propagation (hours-days)
   - Atmospheric drag important
   - Custom force models
   Speed: ★★☆☆☆ (slowest)
   Accuracy: ★★★★★ (highest)
   Duration: Hours-days

DECISION TREE:
==============
1. How long is propagation?
   - < 1 hour: Use Numerical (highest accuracy)
   - 1-5 days: Use DSST or Numerical
   - 1-4 weeks: Use DSST
   - > 1 month: Use DSST or Brouwer-Lyddane
   - Quick estimate: Use Keplerian

2. Is drag significant?
   - Yes: Use Numerical (only option with drag)
   - No: Continue to step 3

3. Do you need full gravity model?
   - Yes (tesseral): Use DSST or Numerical
   - J2-J5 sufficient: Use Brouwer-Lyddane
   - No perturbations: Use Keplerian

4. What's your speed/accuracy requirement?
   - Need 2-10x speedup over numerical: Use DSST
   - Need absolute best accuracy: Use Numerical
   - Need very fast: Use Keplerian or Brouwer-Lyddane
""")

