import numpy as np
import pandas as pd
import sys

# read planetary data
data = pd.read_csv('fake_systems/3_30_0.05.csv', header=None, names=['name', 'p1', 'p2', 'p3', 'p4', 'm1', 'm2', 'm3', 'm4'])

try:
    names = [str(sys.argv[1])]
except:
    names = ["0"]

# for each system
for i in range(len(names)):

    samplename = int(names[i])
    print('system %d' % samplename)
    # finding and isolating the sample data
    # sampleindex=np.where(data["name"]==samplename)[0]
    sampleindex = int(names[i])
    data = data.loc[[sampleindex]]

    Np = 4

    Plist = [data["p1"], data["p2"], data["p3"], data["p4"]]
    Plist = np.array([float(i) for i in Plist])
    mlist = [data["m1"], data["m2"], data["m3"], data["m4"]]
    mlist = np.array([float(i) for i in mlist])

    # add manual masses or periods here (in units of days and earth masses)
    # Plist=np.multiply(1,np.array( [1, 1.82, 3.5, 6.7]))
    # mlist=np.multiply(1,np.array([15,15,15,15]))

    # for resonance checking

    # period ratios of planets
    ratios = np.zeros((Np, Np))
    for x in range(Np):
        for y in range(Np):
            ratios[x, y] = Plist[x] / Plist[y]
    print('normalized periods')
    print(Plist)
    print('period ratios')
    print(ratios)

    # resonances to check for
    orders = 3
    far = 10  # ends at (far+order)/far
    resonances = np.zeros(far * orders)
    resonancesfull = np.zeros((orders, far))
    for x in range(orders):
        for y in range(far):
            resonances[far * x + y] = (y + 2 + x) / (y + 1)
            resonancesfull[x, y] = (y + 2 + x) / (y + 1)

    # print(resonancesfull)
    resonances = np.unique(resonances)

    #check for closest period ratio to resonance
    close2resonance = np.zeros((Np, Np))
    closest2resonance = 10
    for x in range(Np):
        for y in range(Np):
            if y < x:
                value = min(abs(ratios[x, y] - resonances))
                close2resonance[x, y] = value
                closest2resonance = min(value, closest2resonance)

    # tell user results
    if closest2resonance < 0.03:
        print("has resonance, %f \n" % closest2resonance)
    else:
        print("no resonance, %f \n" % closest2resonance)

    # for mutual hill radii separation checking

    # get the star mass
    # Ms=sample["st_mass"][0:1].astype(np.float)
    Ms = 1  # DELETE THIS OR COMMENT OUT IF USING REAL SYSTEMS

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
            if y  <x:
                mut_hill[x, y] = np.mean([a[x], a[y]]) * ((mlist[x] + mlist[y]) * earth / Ms / 3.) ** (1. / 3.)
                separation[x, y] = abs(a[x] - a[y]) / mut_hill[x, y]
                if y == x - 1:
                    singlesep[y] = separation[x, y]
                if y == x - 2:
                    doublesep[y] = separation[x, y]
    maxsinglesep = np.max(singlesep)
    maxdoublesep = np.max(doublesep)

    print('masses')
    print(mlist)
    print('hill radii separations')
    print(separation)
    # assign categories based on separations and report results to user
    # if maxsinglesep>30:
    #   print("too separated, %f"%maxsinglesep)
    # elif maxdoublesep>30:
    #   print("semi-closely packed, %f | %f"%(maxsinglesep, maxdoublesep))
    # else:
    #   print("very closely packed, %f"%maxdoublesep)
