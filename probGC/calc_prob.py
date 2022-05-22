#!/usr/bin/env python
import math


def calc_sigma(method, slope, intercept, r):
    if method == 'devauc':
        lsigma = intercept + (slope*(math.pow(r, 0.25)))
    else:
        lsigma = intercept + (slope*(math.log10(r)))

    sigma = math.pow(10, lsigma)
    return sigma


coeff = [
{'id':1,'name':'NGC4472','dev_a0':3.38,'dev_a1':-1.56,'pow_a0':1.91,'pow_a1':-1.28,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':2,'name':'NGC4406','dev_a0':3.18,'dev_a1':-1.58,'pow_a0':1.66,'pow_a1':-1.24,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':5,'name':'NGC5813','dev_a0':4.18,'dev_a1':-2.23,'pow_a0':2.02,'pow_a1':-1.74,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':6,'name':'NGC3379','dev_a0':2.33,'dev_a1':-1.62,'pow_a0':0.83,'pow_a1':-1.41,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':8,'name':'NGC4594','dev_a0':3.76,'dev_a1':-2.11,'pow_a0':1.87,'pow_a1':-1.85,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':9,'name':'NGC4762','dev_a0':4.23,'dev_a1':-3.09,'pow_a0':1.1,'pow_a1':-1.86,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':10,'name':'NGC1023','dev_a0':3.41,'dev_a1':-2.15,'pow_a0':1.23,'pow_a1':-1.36,'method':'devauc','contam':1.1,'gclf':0.57,'extent':15.0},
{'id':11,'name':'NGC5866','dev_a0':3.52,'dev_a1':-2.52,'pow_a0':0.96,'pow_a1':-1.76,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':12,'name':'NGC7332','dev_a0':3.24,'dev_a1':-2.32,'pow_a0':0.93,'pow_a1':-1.6,'method':'pow','contam':1.1,'gclf':0.36,'extent':2.5},
{'id':13,'name':'NGC4754','dev_a0':3.44,'dev_a1':-2.61,'pow_a0':0.78,'pow_a1':-1.41,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':14,'name':'NGC3384','dev_a0':2.06,'dev_a1':-1.51,'pow_a0':0.65,'pow_a1':-1.29,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':15,'name':'NGC7457','dev_a0':4.09,'dev_a1':-3.01,'pow_a0':1.03,'pow_a1':-1.82,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':17,'name':'NGC7814','dev_a0':2.86,'dev_a1':-1.81,'pow_a0':1.03,'pow_a1':-0.98,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':18,'name':'NGC3556','dev_a0':4.48,'dev_a1':-3.13,'pow_a0':1.43,'pow_a1':-2.27,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':19,'name':'NGC2683','dev_a0':4.28,'dev_a1':-3.39,'pow_a0':0.93,'pow_a1':-2.35,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':20,'name':'NGC4157','dev_a0':2.02,'dev_a1':-1.71,'pow_a0':0.34,'pow_a1':-1.27,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':21,'name':'NGC7331','dev_a0':2.33,'dev_a1':-1.57,'pow_a0':0.81,'pow_a1':-1.16,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':24,'name':'NGC4013','dev_a0':4.02,'dev_a1':-3.11,'pow_a0':0.88,'pow_a1':-1.73,'method':'devauc','contam':1.1,'gclf':0.66,'extent':15.0},
{'id':25,'name':'NGC891','dev_a0':3.02,'dev_a1':-2.53,'pow_a0':0.45,'pow_a1':-1.62,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':26,'name':'NGC7339','dev_a0':2.89,'dev_a1':-1.78,'pow_a0':1.09,'pow_a1':-1.24,'method':'pow','contam':1.1,'gclf':0.33,'extent':1.5},
{'id':27,'name':'NGC1055','dev_a0':3.43,'dev_a1':-2.6,'pow_a0':0.77,'pow_a1':-1.55,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':34,'name':'NGC4649','dev_a0':3.42,'dev_a1':-1.67,'pow_a0':1.64,'pow_a1':-1.03,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':42,'name':'NGC1172','dev_a0':4.58,'dev_a1':-3.1,'pow_a0':1.44,'pow_a1':-1.82,'method':'devauc','contam':1.1,'gclf':0.5,'extent':15.0},
{'id':43,'name':'NGC4382','dev_a0':4.14,'dev_a1':-2.37,'pow_a0':1.73,'pow_a1':-1.57,'method':'pow','contam':1.1,'gclf':0.37,'extent':20.5},
{'id':45,'name':'NGC4621','dev_a0':3.46,'dev_a1':-2.01,'pow_a0':1.41,'pow_a1':-1.35,'method':'devauc','contam':1.1,'gclf':0.36,'extent':14.2},
{'id':54,'name':'NGC5846','dev_a0':3.61,'dev_a1':-1.80,'pow_a0':1.77,'pow_a1':-1.20,'method':'devauc','contam':1.1,'gclf':0.36,'extent':14.2}]

cg = coeff[0]

f = open('gcc.dat', 'r')
for l in f.readlines():
    vals = l.split()
    id = int(vals[0])
    r = float(vals[1])
    galaxy = int(vals[2])

    if galaxy != cg['id']:
        for c in coeff:
            if c['id'] == galaxy:
                cg = c

    if galaxy != cg['id']:
        continue 
    if cg['method'] == 'devauc':
        intercept = cg['dev_a0']
        slope = cg['dev_a1']
    else:
        intercept = cg['pow_a0']
        slope = cg['pow_a1']

    sigma = calc_sigma(cg['method'], slope, intercept, r)

    corr_sigma = sigma * cg['gclf']

    prob = corr_sigma / (corr_sigma + cg['contam'])

    if prob < 0:
        prob = 0.0

    print r/cg['extent'], prob, cg['name']

