import glob
import sys
import datetime
from os import path

narg =  len(sys.argv)
folder = "./data"
ext = "log"
ntot = 50
print(narg)
if 1 < narg:
    folder = sys.argv[1]
    if 2 < narg:
        ext = sys.argv[2]
        if 3 < narg:
            ntot = int(sys.argv[3])

files = [f for f in glob.glob(folder + "**/*." + ext, recursive=False)]

# Load all log files
print("Loading data...")
min_time = 1E10
min_time = 1E10
users = {}
count = 0
user_index = {}
user_outcome = [0] * ntot
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
    
    # Events are stored last to first in log files, reverting the order
    events.reverse()
    users[case_id] = events
    user_index[case_id] = count
    count += 1

# Construct infection network
ninf = 0
nknown = 0
nmiss = 0
nsurv = 0
ndead = 0
npeer = 0
inf_network = []
for key in users:
#     print("=================", key)
    events = users[key]
    if len(events) == 0: 
        continue
#         ninf += 1
    pkey = ""
    infect = None
    has_inf_event = False
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
                user_outcome[user_index[key]] = 1
                nsurv += 1
            elif data == "DEAD":
                user_outcome[user_index[key]] = 2                
                ndead += 1
            if not infect:
                # Infection edge without origin
                infect = [{"a":"unk", "b":key, "t":date, "s":"-"}]
            ninf += 1
        elif ev["type"] == "INF":
            if pkey: 
                # Duplicated parent, skipping
                continue
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
                    infect = [{"a":pkey, "b":key, "t":date, "s": strain}]
                    npeer += 1
                    has_inf_event = True
                else:
                    infect = [{"a":"unk", "b":key, "t":date, "s": strain}]
                    has_inf_event = False                 
            elif "CASE0" in data:
                strain = data[6:-1]
                # print(strain, "0", "->", key)
                infect = [{"a":"zero", "b":key, "t":date, "s": strain}]
                has_inf_event = True
            elif "SOURCE" in data:
                infect = [{"a":"src", "b":key, "t":date, "s": strain}]
                has_inf_event = True

    if has_inf_event:
        nknown += 1
    else:
        nmiss += 1
        
    inf_network += infect

print("Total number of users:", ntot)
print("Total number of cases:", ninf)
print("Total number of deaths:", ndead)
print("Total number of survivors:", nsurv)
print("Number of infections with known source:", nknown)
print("Number of infections from peer:", npeer)
print("Number of infections with missing source:", nmiss)

# Plotting

from igraph import *

g = Graph(directed=True)
g.add_vertices(ntot)

g.vs["outcome"] = user_outcome
color_dict = {0: "Green Yellow", 1: "Deep Sky Blue", 2: "Orange Red"}
g.vs["color"] = [color_dict[out] for out in g.vs["outcome"]]

for edge in inf_network:
    n0 = edge["a"]
    n1 = edge["b"]
    if n0 in user_index and n1 in user_index:
        print(user_index[n0], "->", user_index[n1])
        g.add_edges([(user_index[n0], user_index[n1])])
        
print(g)

layout = g.layout("fr")
plot(g, layout = layout, vertex_size = 7, edge_arrow_width = 1, edge_arrow_size = 0.5)
