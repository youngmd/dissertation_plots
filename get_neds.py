from astroquery.ned import Ned
import urllib2
from BeautifulSoup import BeautifulSoup
import numpy as np
import sys
import os
import math
import csv

infile = sys.argv[1]

# f = open(infile, 'r')
# import csv


Ks = []
Kerrs = []

Vs = []

with open(infile, 'rb') as f:
    reader = csv.DictReader(f)
    for mrow in reader:
        try:
            name = mrow['Galaxy']


# for t in targets:
#     try:
#
#         name = t.strip()
        # name = v[0] + " " + v[1]
        # mass = float(v[2])
        # mass_err = float(v[3])
        # mass_str = "%.1f$\\pm$%.1f" % (mass, mass_err)
        # # mass_str = "%.1f" % (mass)
        #
        # density = float(v[4])
        # re = float(v[5])
        # re_err = float(v[6])
        # re_str = "%.1f$\\pm$%.1f" % (re,re_err)
        #
        # sige = v[7]
        # sige_err = v[8]
        # sige_str = "%s$\\pm$%s" % (sige, sige_err)
        #
        # lgmbh = float(v[9])
        # lgmbhup = float(v[10])
        # lgmbhdown = float(v[11])
        # lgstr = "$%.2f_{%.2f}^{%.2f}$" % (lgmbh, lgmbhup, lgmbhdown)
        #
        # if lgstr.startswith('$999'):
        #     lgstr = ''

            flags = []
            # Get RC3 classification from NED
            ut = name.replace(" ","")
            url = "https://ned.ipac.caltech.edu/cgi-bin/NEDatt?objname="+ut
            # print url
            html = urllib2.urlopen(url).read()
            bs = BeautifulSoup(html)
            tr = bs.find(lambda tag: tag.name=='tr' and 'Galaxy Morphology' in tag.text.encode('utf-8')).findNext('tr')
            cells = tr.findChildren('td')
            morph = cells[0].text
            morph = morph.replace('edge-on', '')
            morph = morph.translate({ord(i):None for i in '()+:?rsR'})
            morph = morph.replace('-', ' ')
            morph = morph.replace('00', '0')
            morph = morph.replace(' ', '')
            morph = morph.strip()

            # Get Distance from NED
            tr = bs.find(lambda tag: tag.name=='tr' and 'SBF' in tag.text.encode('utf-8') and '2001ApJ...546..681T' in tag.text.encode('utf-8'))
            dflag = ''

            if not tr:
                tr = bs.find(lambda tag: tag.name=='tr' and 'Tully-Fisher' in tag.text.encode('utf-8') and '2009ApJS' in tag.text.encode('utf-8'))
                dflag = 'a'

            if not tr:
                tr = bs.find(lambda tag: tag.name=='tr' and 'GCLF' in tag.text.encode('utf-8') and '2007ApJS..171..101J' in tag.text.encode('utf-8'))
                dflag = 'b'

            if not tr:
                tr = bs.find(lambda tag: tag.name == 'tr' and 'Tully' in tag.text.encode('utf-8') and '465...71T' in tag.text.encode('utf-8'))
                dflag = 'c'

            if not tr:
                tr = bs.find(lambda tag: tag.name == 'tr' and 'Tully' in tag.text.encode('utf-8') and '2009AJ....138..323T' in tag.text.encode('utf-8'))
                dflag = 'd'

            if not tr:
                tr = bs.find(lambda tag: tag.name == 'tr' and 'Tully' in tag.text.encode('utf-8') and '1988NBGC' in tag.text.encode('utf-8'))
                dflag = 'e'

            if not tr:
                dmod_str = 'NA'
                dist = 'NA'
            else:
                cells = tr.findChildren('td')
                distance = cells[6].text
                d = distance.split()
                dmod = float(d[2])
                dmod_err = float(d[4])
                dist = float(d[7])
                dmod_str = "%.2f$\\pm$%.2f" % (dmod, dmod_err)

            if dflag:
                flags.append(dflag)

            # Get photometry from NED
            rt = Ned.get_table(ut, table="photometry")
            Vflag = ''
            rowV  = np.where(rt['Observed Passband'] == 'V (V_T^0)')
            rowV2 = np.where(rt['Observed Passband'] == 'V (V_T)')
            if not rowV[0].size or not rowV2[0].size:
                Vstr = 'NA'
                Vflag = 'f'
                Vs.append('NA')
            else:
                V = rt[rowV]['Photometry Measurement'].data.data[0]
                V_err = float(rt[rowV2]['Uncertainty'].data.data[0][3:])
                V = V - dmod
                Vs.append(V)
                V_err = math.sqrt(math.pow(V_err, 2) + math.pow(dmod_err, 2))
                Vstr = "%.2f$\\pm$%.2f" % (V, V_err)
                # Vstr = "%.1f" % (V)
            if Vflag:
                flags.append(Vflag)

            Kflag = ''
            row = np.where(rt['Observed Passband'] == 'K_tot (2MASS LGA)')
            row2 = np.where(rt['Observed Passband'] == 'K_s_total (2MASS)')
            if not row[0].size and not row2[0].size:
                    # print "Could not find 2mass measurements for %s" % t
                    K = 'NA'
                    K_err = 'NA'
            elif not row[0].size:
                row = row2
                K = rt[row]['Photometry Measurement'].data.data[0]
                K_err = float(rt[row]['Uncertainty'].data.data[0][3:])
            else:
                K = rt[row]['Photometry Measurement'].data.data[0]
                K_err = float(rt[row]['Uncertainty'].data.data[0][3:])
            if K == 'NA':
                Kstr = 'NA'
                Kflag = 'g'
            else:
                K = K - dmod
                K_err = math.sqrt(math.pow(K_err, 2) + math.pow(dmod_err, 2))
                Kstr =  "%.2f$\\pm$%.2f" % (K, K_err)
                # Kstr =  "%.1f" % (K)
            if Kflag:
                flags.append(Kflag)

            if flags:
                flags = '$,$'.join(flags)
                flags = '$' + flags + '$'
            else:
                flags = ''

            print name,mrow["dmod"],dmod_str,morph,Vstr,Kstr,flags
        # print "%s & %s & %s & %s & %s & %s & %.2f & %s & %s & %s & %s \\\\" % (name,
        #                                         morph, mass_str, Vstr,
        #                                         Kstr, dmod_str, density,
        #                                         re_str, sige_str, lgstr,
        #                                         flags)
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)
            print e
            print "%s - NA" % t

for V in Vs:
    print V

# print "======="

# for K in Kerrs:
#     print K