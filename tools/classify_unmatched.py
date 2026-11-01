"""Exploratory, outside the held-out protocol (written after the Greek test,
2026-11-01): where the unmatched alarms of one ayzek report locate.

    python3 tools/classify_unmatched.py REPORT LAT LON

REPORT is an ayzek report (`--record` run output), LAT LON the first event of the
sequence.
A location counts if it is a P-time location from >=5 stations with a 68% radius
<=20 km, or a pick location from >=3 stations with rms <=1 s (the last of each).
Classes: 'sequence' (counted location within 30 km of the sequence's first event),
'elsewhere' (counted location farther), 'unlocated' (no counted location)."""
import re, sys, math, collections
rep, lat0, lon0 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
t = open(rep).read()
def km(a,b,c,d):
    p=math.radians; return 6371*2*math.asin(math.sqrt(math.sin(p(c-a)/2)**2+math.cos(p(a))*math.cos(p(c))*math.sin(p(d-b)/2)**2))
alarms = {}
for m in re.finditer(r'^ALARM\s+#(\d+)\s+(\S+)\s+earthquake detected by ([^(\n]+?)\s*(\(AFAD.*\))?$', t, re.M):
    alarms[int(m.group(1))] = dict(time=m.group(2), sta=[x.strip() for x in m.group(3).split(',')], matched=bool(m.group(4)), loc=None)
for m in re.finditer(r'^PLOC\s+#(\d+).*locate it at ([\d.]+)N ([\d.]+)E \+-(\d+) km.*from (\d+) stations', t, re.M):
    i=int(m.group(1)); 
    if int(m.group(5))>=5 and int(m.group(4))<=20: alarms[i]['loc']=(float(m.group(2)),float(m.group(3)),'P')
for m in re.finditer(r'^LOCATE\s+#(\d+).*located at ([\d.]+)N ([\d.]+)E.*rms ([\d.]+) s from (\d+) stations', t, re.M):
    i=int(m.group(1))
    if int(m.group(5))>=3 and float(m.group(4))<=1.0: alarms[i]['loc']=(float(m.group(2)),float(m.group(3)),'picks')
# unmatched list from the summary
um = set(int(x) for x in re.findall(r'^\s+\S+ UTC\s+#(\d+)\s', t.split('False alarms')[1].split('---')[0], re.M)) if 'False alarms' in t else set()
cls = collections.Counter(); near=collections.Counter(); sta=collections.Counter()
for i in um:
    a=alarms[i]
    for s in a['sta']: sta[s]+=1
    if a['loc'] is None: cls['unlocated']+=1
    else:
        d=km(lat0,lon0,a['loc'][0],a['loc'][1]); c='sequence' if d<=30 else 'elsewhere'; cls[c]+=1; near[(c,a['loc'][2])]+=1
print(f"alarms {len(alarms)}, matched {sum(a['matched'] for a in alarms.values())}, unmatched {len(um)}: "+", ".join(f"{k} {v}" for k,v in cls.most_common())+f"  | by method {dict(near)}")
print("  stations in unmatched alarms:", ", ".join(f"{k} {v}" for k,v in sta.most_common(8)))
