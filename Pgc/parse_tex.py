import sys

infile = sys.argv[1]

f = open(infile, 'r')
lines = f.readlines()
for l in lines:
    vals = l.split('&')
    vmag = vals[3]
    r = vals[8]
    v = vals[10]
    yn = vals[11]
    v = v.split('$\\pm$')
    v_r = v[0].strip()
    if v_r.startswith('-'):
        v_r = v_r[1:]
    v_err = v[1].strip()
    vmag = vmag.strip()
    r = r.strip()
    yn = yn[:2]
    yn = yn.strip()
    print vmag,r,v_r,v_err,yn
