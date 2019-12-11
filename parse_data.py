import glob
import sys
import datetime
from os import path

narg =  len(sys.argv)
folder = "./data"
ext = "log"
if 1 < narg:
    folder = sys.argv[1]
    if 2 < narg:
        ext = sys.argv[2]

files = [f for f in glob.glob(folder + "**/*." + ext, recursive=False)]

print("Loading data...")
min_time = 1E10
users = {}
for fn in files:
    bname = path.basename(fn)
    idx = bname.rfind('-')
    case_id = bname[0:idx]
    events = []
    with open(fn, 'r') as f:
        lst = list(enumerate(f))
        n = len(lst)

        for i, line in lst:
            if i == 0 or i == n - 1:
                continue
            line = line.strip()
            line = line[1:-2]
            parts = line.split(',')
            time = int(parts[0].split(':')[1])

            evstr = parts[1]
            idx = evstr.find(':')
            evtyp = evstr[0:idx]
            if idx < len(evstr):
                evdat = evstr[idx+1:]
            else:
                evdat = None

            # print(time, evtyp, evdat)
            events += [{"time": time, "type": evtyp, "data": evdat}]

            min_time = min(min_time, time)
    
    events.reverse()
    users[case_id] = events

# Normalize time
ntot = len(users)
ninf = 0
nmiss = 0
nsurv = 0
ndead = 0
inf_network = []
for key in users:
    print("=================", key)
    events = users[key]
    if 0 < len(events):
        ninf += 1
    for ev in events:
        # ev["time"] -= min_time
        # ev["time"] = ev["time"] - min_time
        date = datetime.datetime.fromtimestamp(ev["time"])
        # date = date.strftime('%Y-%m-%d %H:%M:%S')
        date = date.strftime('%H:%M:%S')
        # print(ev["time"] - min_time, date, ev["type"], ev["data"])
        data = ev["data"]
        if ev["type"] == "OUT":
            if data == "RECOVERED":
                nsurv += 1
            elif data == "DEAD":
                ndead += 1
        elif ev["type"] == "INF":
            if "PEER" in data:
                pstr = data[5:-1]
                if ":" in pstr:
                    pieces = pstr.split(":")
                    pkey = pieces[0]
                    strain = pieces[1]
                else:
                    pkey = pstr
                    strain = "0"
                if pkey in users:
                    # print(strain, pkey, "->", key)
                    inf_network += [{"a":pkey, "b":key, "t":date, "s": strain}]
                else:
                    nmiss += 1
            elif "CASE0" in data:
                strain = data[6:-1]
                # print(strain, "0", "->", key)
                inf_network += [{"a":"0", "b":key, "t":date, "s": strain}]
            elif "SOURCE" in data:
                continue

print("Total number of users:", ntot)
print("Total number of cases:", ninf)
print("Total number of deaths:", ndead)
print("Total number of survivors:", nsurv)
print("Number of infections with known source:", len(inf_network))
print("Number of infections with missing source:", nmiss)

