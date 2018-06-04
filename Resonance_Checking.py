import numpy as np
import pandas as pd
import make_4planet as m4

# read planetary data
data = pd.read_csv('planets_mod.csv', header=40)

names = list(set(data["pl_hostname"]))  # list of unique names
names.sort()  # sort the list

# do you want to check for resonances or for mutual hill radii separation?
# NEEDS TO BE CHANGED EVERYTIME!
checkresonance = 1

names = ["Ari Fake 10 0.1"]

# for each system
for i in range(len(names)):

    samplename = names[i]
    print(samplename)

    # finding and isolating the sample data
    sampleindex = np.where(data["pl_hostname"] == samplename)[0]
    sample = data[sampleindex[0]:sampleindex[-1] + 1]
    # sort the samples in order of their periods
    sample = sample.sort_values(by="pl_orbper")

    # number of planets
    Np = sample.shape[0]

    Plist = np.zeros(Np)
    mlist = np.zeros(Np)

    # how many samples to check over (generated values are distributions)
    n_sims = 100

    if samplename not in ["HR 8799", "K2-72", "Kepler-176", "Kepler-341", "Kepler-37"]:
        for i in range(Np):
            mass, period, abc1, abc2, abc3 = m4.planetparameters(sample, i, n_sims)

            # maybe median should be used?
            mlist[i] = np.mean(mass)
            Plist[i] = np.mean(period)

    # add manual masses or periods here (in units of days and earth masses)
    # Plist=np.multiply(1,np.array( [1, 1.82, 3.5, 6.7]))
    # mlist=np.multiply(1,np.array([23,23,23,23]))

    # for resonance checking
    if checkresonance:

        # period ratios of planets
        ratios = np.zeros((Np, Np))
        for x in range(Np):
            for y in range(Np):
                ratios[x, y] = Plist[x] / Plist[y]
        print(ratios)

        # resonances to check for
        orders = 5
        far = 10  # ends at (far+order)/far
        resonances = np.zeros(far * orders)
        resonancesfull = np.zeros((orders, far))
        for x in range(orders):
            for y in range(far):
                resonances[far * x + y] = (y + 2 + x) / (y + 1)
                resonancesfull[x, y] = (y + 2 + x) / (y + 1)

        print(resonancesfull)
        resonances = np.unique(resonances)

        # check for closest period ratio to resonance
        close2resonance = np.zeros((Np, Np))
        closest2resonance = 10
        for x in range(Np):
            for y in range(Np):
                if y < x:
                    value = min(abs(ratios[x, y] - resonances))
                    close2resonance[x, y] = value
                    closest2resonance = min(value, closest2resonance)

        # tell user results
        if closest2resonance < 0.05:
            print("has resonance, %f" % closest2resonance)
        else:
            print("no resonance, %f" % closest2resonance)

    # for mutual hill radii separation checking
    else:

        # get the star mass
        Ms = sample["st_mass"][0:1].astype(np.float)
        # Ms = 1 #DELETE THIS OR COMMENT OUT IF USING REAL SYSTEMS

        # get semi-major axes
        earth = 0.000003003
        a = np.zeros(Np)
        for i in range(Np):
            a[i] = ((Plist[i] / 365 * (2 * np.pi)) ** 2 * Ms) ** (1. / 3.)

        # find mutual hill radii separations and record those that are single
        # separated and doubly separated (one planet between the pair)
        # (this could be generalized)
        mut_hill = np.zeros((Np, Np))
        separation = np.zeros((Np, Np))
        singlesep = np.zeros(Np - 1)
        doublesep = np.zeros(Np - 2)
        for x in range(Np):
            for y in range(Np):
                if y < x:
                    mut_hill[x, y] = np.mean([a[x], a[y]]) * ((mlist[x] + mlist[y]) * earth / Ms / 3.) ** (1. / 3.)
                    separation[x, y] = abs(a[x] - a[y]) / mut_hill[x, y]
                    if y == x - 1:
                        singlesep[y] = separation[x, y]
                    if y == x - 2:
                        doublesep[y] = separation[x, y]
        maxsinglesep = np.max(singlesep)
        maxdoublesep = np.max(doublesep)

        print(separation)
        # assign categories based on separations and report results to user
        if maxsinglesep > 30:
            print("too separated, %f" % maxsinglesep)
        elif maxdoublesep > 30:
            print("semi-closely packed, %f | %f" % (maxsinglesep, maxdoublesep))
        else:
            print("very closely packed, %f" % maxdoublesep)
        print(mlist)
