import numpy as np
from scipy import optimize


# This function returns the combined mass for two planets
# needed to be at a certain hill radius limit.
# r is the period ratio between the two planets
# mass_limit scales the result. Effectively used to
# subtract the value of one planet to find a limit on the other
def f(r, hill_limit, mass_limit):
    Mstar = Ms / earth
    return ((24 * Mstar * (r ** (2 / 3) - 1) ** 3) / (hill_limit ** 3 * (r ** (2 / 3) + 1) ** 3) - mass_limit)


# Draws a valid period that is not close to a resonance and within the passed limits
def next_period(p, index, lower_ratio_limit, upper_ratio_limit):
    closest2resonance = 0
    j = index + 1

    # if the current period is too close to a resonance, draw another
    while closest2resonance < resonance_tolerance:

        # draw the period within the given limits
        p[j - 1] = p[j - 2] * 10 ** (np.random.uniform(low=np.log10(lower_ratio_limit), high=np.log10(upper_ratio_limit)))

        # period ratios of planets
        ratios = np.zeros((j, j))
        for x in range(j):
            for y in range(j):
                ratios[x, y] = p[x] / p[y]

        # check for closest period ratio to resonance
        close2resonance = np.zeros((j, j))
        for x in range(j):
            for y in range(j):
                if y < x:

                    #if you are close in, only worry about
                    # first order resonances
                    if ratios[x, y] < 1.5:
                        close2resonance[x, y] = min(abs(ratios[x, y] - first_order_resonances))
                    else:
                        close2resonance[x, y] = min(abs(ratios[x, y] - resonances))

        # find the smallest nonzero value in the close to resonance matrix
        closest2resonance = min(i for i in close2resonance.flatten() if i > 0)

    # return the period ratio (for convenient printing later)
    return ratios

# only seems to work for 1 call. otherwise the optimize function freaks out
# Traceback (most recent call last):
#   File "generate_compact_system.py", line 142, in <module>
#     lower_ratio_limit = max(min_resonance_period_ratio,optimize.brentq(f, 1, 10, args=(lowerhill,upperplanetmass+m[i-1])))
#   File "C:\Users\chris\Anaconda3\lib\site-packages\scipy\optimize\zeros.py", line 519, in brentq
#     r = _zeros._brentq(f,a,b,xtol,rtol,maxiter,args,full_output,disp)
# TypeError: '_io.TextIOWrapper' object is not callable
try:
    # repeat=int(sys.argv[1])
    repeat = 1
except:
    repeat = 1

# Setting some constants
earth = 0.000003003  # Earth mass/Sun mass
Ms = 1  # Mass of star for the test system in solar masses

# Number of planets we want in the system
Np = 4

# Setting hard limits
lowerplanetmass = 1 / 30  # Earth masses
lowerplanetmass = 2 / 3  # Earth masses
# lowerplanetmass=32/8 #Earth masses (1/8 of minimum system mass)
# bullymass=12 #Earth masses

# upperplanetmass=30  # Earth masses
upperplanetmass = 3 / 2  # Earth masses

lowerhill = 1  # mutual hill radii
upperhill = 15  # mutual hill radii
secondaryhill = 30  # mutual hill radii
resonance_tolerance = 0.01  # how close to a resonance a period ratio can be
maxpr = 1.5  # maximum period ratio separation

# resonances to check for
orders = 3
far = 10  # ends at (far+order)/far
resonances = np.zeros(far * orders + 1)
# resonances[-1]=1 #added for printing purposes?
for x in range(orders):
    for y in range(far):
        resonances[far * x + y] = (y + 2 + x) / (y + 1)
resonances = np.unique(resonances)
# resonances=np.delete(resonances, 0) #remove 1 from the list
# print(resonances)

# getting only the first order resonances
first_order_resonances = np.array([(2 + x) / (1 + x) for x in range(far)])

# find the smallest period ratio allowed by the closeness to resonance
# criterion
i = 0
min_resonance_period_ratio = first_order_resonances[0]
while i < len(first_order_resonances):

    # If the current minimum period ratio is too close to the next highest
    # resonance, increase it
    if min_resonance_period_ratio >= first_order_resonances[i] - resonance_tolerance:
        min_resonance_period_ratio = first_order_resonances[i] + resonance_tolerance
        i += 1

    # if it is far enough from the next highest resonance, leave it where it is
    # and end the loop
    else:

        # break the loop
        i = len(first_order_resonances)
# print(min_resonance_period_ratio)

for asldnlkasjd in range(repeat):

    # Initializing period and mass lists
    P = np.zeros(Np)
    P[0] = 1
    m = np.zeros(Np)

    # draw periods and masses iteratively until you get a working result
    i = 0
    resets = 0
    while i < Np:

        # for the first draw, you just need to get a mass within the hard
        # limits
        if i == 0:

            m[i] = 10 ** np.random.uniform(low=np.log10(lowerplanetmass), high=np.log10(upperplanetmass))
            i += 1

        # for other draws, get period and mass
        elif i > 0:

            # get limits on period ratios from constraints on total mass and period ratios
            lower_ratio_limit = max(min_resonance_period_ratio, optimize.brentq(f, 1, maxpr, args=(lowerhill, upperplanetmass + m[i - 1])))
            upper_ratio_limit = optimize.brenth(f, 1, maxpr, args=(upperhill, upperplanetmass + m[i - 1]))

            # if the constraints are impossible, start fresh
            if lower_ratio_limit > upper_ratio_limit:
                i = 0
                resets += 1
                print("reset b/c of period constraints")

            # if a period can be drawn, draw it and get a mass
            else:
                ratios = next_period(P, i, lower_ratio_limit, upper_ratio_limit)

                #after the second planet, look at secondary hill radii separations as well
                if i > 1:
                    lower_mass_limit = max(lowerplanetmass, f(P[i] / P[i - 1], upperhill, m[i - 1]), f(P[i] / P[i - 2], secondaryhill, m[i - 2]))

                else:
                    lower_mass_limit = max(lowerplanetmass, f(P[i] / P[i - 1], upperhill, m[i - 1]))
                upper_mass_limit = min(upperplanetmass, f(P[i] / P[i - 1], lowerhill, m[i - 1]))

                # if i==3:
                #   lower_mass_limit=max(bullymass, f(P[i]/P[i-1],upperhill,m[i-1]), f(P[i]/P[i-2],secondaryhill,m[i-2]))
                #   upper_mass_limit=min(30, f(P[i]/P[i-1],lowerhill,m[i-1]))

                # if the constraints are impossible, start fresh
                if lower_mass_limit > upper_mass_limit:
                    i = 0
                    resets += 1
                    print("reset b/c of mass constraints")
                # if a mass can be drawn, draw it
                else:
                    m[i] = 10 ** np.random.uniform(low=np.log10(lower_mass_limit), high=np.log10(upper_mass_limit))

                    # move on to the next planet
                    i += 1

    # the smallest orbital period is 1/10 of a year (in days)
    scalefactor = 365.25 / 10
    Preal = P * scalefactor  # I'm pretty sure this doesn't matter

    # get semi-major axes
    a = np.zeros(Np)
    for i in range(Np):
        a[i] = ((Preal[i] / 365 * (2 * np.pi)) ** 2 * Ms) ** (1. / 3.)

    # find mutual hill radii separations for printing
    mut_hill = np.zeros((Np, Np))
    separation = np.zeros((Np, Np))
    for x in range(Np):
        for y in range(Np):
            if y < x:
                mut_hill[x, y] = np.mean([a[x], a[y]]) * ((m[x] + m[y]) * earth / Ms / 3.) ** (1. / 3.)
                separation[x, y] = abs(a[x] - a[y]) / mut_hill[x, y]

    print("system made")
    print("and it took %d resets" % (resets))
    # print("\nPeriods (scaled to first period)\n")
    # print(P)
    # print("\nPeriod Ratios\n")
    # print(ratios)
    # print("\nPlanet Masses (Earth masses)\n")
    # print(m)
    # print("\nHill Separations\n")
    # print(separation)

    csvname = 'fake_systems/%1.1f_%1.1f_%1.2f.csv' % (lowerplanetmass, upperplanetmass, resonance_tolerance)
    # csvname='fake_systems/%d_%d_%1.2f_small.csv'%(lowerhill,upperhill,resonance_tolerance)
    try:
        fo = open(csvname, "r")
        lines = fo.readlines()[-1]
        index = int(lines.split(',')[0]) + 1
        fo.close()
    except:
        index = 0

    f = open(csvname, "a")
    f.write('%d, %e, %e, %e, %e, %e, %e, %e, %e \n' % (index, P[0], P[1], P[2], P[3], m[0], m[1], m[2], m[3]))
    f.close()
