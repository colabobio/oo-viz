import os, sys, json
from os import path
from datetime import datetime, timedelta, date
import pytz

from igraph import *

from PIL import Image, ImageDraw, ImageFont

import seaborn as sns
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as clr

# Load properties
if len(sys.argv) < 1:
    print("JSON files with simulation properties are missing")
    exit(1)

json_fname = sys.argv[1]
with open(json_fname) as f:
    props = json.load(f)

title = props["title"]
base_folder = props["base_folder"]
sim_id = props["sim_id"]
sim_tz = props["sim_tz"]
time_step_min = props["time_step_min"]

if "time0" in props and "time1" in props:
    time0 = props["time0"]
    time1 = props["time1"]
else:
    time0 = time1 = ''

if "use_new_id_schema" in props:
    use_new_id_schema = props["use_new_id_schema"]
else:
    use_new_id_schema = False

# Configuration

# Coded status:
# https://matplotlib.org/3.1.0/gallery/color/named_colors.html
status_color = {0: (1, 1, 1, 0),              # Susceptible (transparent)
                1: clr.to_hex("darkorange"),  # Infected (index case)
                2: clr.to_hex("darkorange"),  # Infected (from someone else)
                3: clr.to_hex("darkgrey"),        # Dead 
                4: clr.to_hex("mediumseagreen"),  # Recovered 
                5: clr.to_hex("darkorchid")       # Vaccinated                
               } 

# https://github.com/google/fonts/tree/master/apache
label_font = ImageFont.truetype("Roboto-Regular.ttf", size=24)

# Visual style of the infection network
# https://igraph.org/python/versions/latest/tutorial.html#vertex-attributes-controlling-graph-plots
style = {}
style["bbox"] = (1200, 800)
style["margin"] = 15
style["vertex_size"] = 7
style["vertex_frame_width"] = 0
style["vertex_label_size"] = 5
style["edge_arrow_size"] = 0.6
style["edge_arrow_width"] = 0.6
style["edge_curved"] = False

data_folder = path.join(base_folder, "data")
output_root = path.join(base_folder, "output")
output_folder = path.join(output_root, "infections")
movie_folder = path.join(output_root, "movies")
if not path.exists(output_folder):
    os.makedirs(output_folder)
if not path.exists(movie_folder):
    os.makedirs(movie_folder)

# Print warning messages to the console when parsing data
print_data_warnings = True
    
frame_format = "png"

# Time delta for plots in seconds
time_delta_sec = 60 * time_step_min

# Parameters of the layout algorithm, the anim_steps is how many times the fruchterman-reingold (fr)
# algorithm is run per time delta, the higher the smoother the animation will be.
# fr_niter controls the number of iterations to perform by the fr algorithm.
# The product of these two numbers should be around 200 ~ 500
anim_steps_per_time_delta = 30
fr_niter = 10

# https://howchoo.com/g/ywi5m2vkodk/working-with-datetime-objects-and-timezones-in-python
# https://itnext.io/working-with-timezone-and-python-using-pytz-library-4931e61e5152
timezone = pytz.timezone(sim_tz)

if time0 and time1:
    obs_date0 = timezone.localize(datetime.strptime(time0, '%b %d %Y %I:%M%p'))
    obs_date1 = timezone.localize(datetime.strptime(time1, '%b %d %Y %I:%M%p'))
else:
    obs_date0 = None
    obs_date1 = None

    # Some utility functions

def get_infection_list(events):
    infections = events[(events["type"] == "infection")]
    
    ilist = []
    itimes = {}
    infected = infections.user_id.values
    peers = infections.inf.values
    timestamp = infections.time.values
    for id1, peer0, ts in zip(infected, peers, timestamp):
        n1 = user_index[id1]
            
        if "PEER" in peer0:
            if use_new_id_schema:
                # New schema
                id0 = int(peer0[peer0.index("[") + 1:peer0.index(":")])
                if id0 in user_index:
                    n0 = user_index[id0]                    
                    add_infection = True
                    for e in ilist:
                        if e[1] == n1:
                            pid0 = index_user[e[0]]
                            ts0 = itimes[(pid0, id1)]
                            if abs(ts - ts0) <= time_delta_sec:
                                add_infection = False
                                if pid0 == id0:                                
                                    print("Duplicated infection:", id1, "was already infected by", id0, "in the last", time_step_min, "minutes")
                                else:
                                    print("Multiple infection:", id1, "is being infected by", id0, "but was already infected by", pid0, "in the last", time_step_min, "minutes")
                                break    

                    if add_infection: 
                        ilist += [(n0, n1)]
                        itimes[(id0, id1)] = ts
                elif print_data_warnings:
                    print("Cannot find peer", id0)                    
            else:    
                # Old schema (sims before 2022): p2p id is in the infection column
                p2p0 = peer0[peer0.index("[") + 1:peer0.index(":")]
                if p2p0 in p2pToId:
                    id0 = p2pToId[p2p0]
                    if id0 in user_index:
                        n0 = user_index[id0]
                        if not (n0, n1) in ilist:                        
                            ilist += [(n0, n1)]
                        elif print_data_warnings:
                            print("Duplicated infection", id0, id1)  
                elif print_data_warnings:
                    print("Cannot find peer", p2p0)                        
            
    return ilist    

def get_node_status(events, status0 = None):
    if status0 == None:
         status = [0] * len(users)
    else:            
        status = status0

    inf = events[events["type"] == "infection"]
    infMap = pd.Series(inf.inf.values, index=inf.user_id).to_dict()
    for kid in infMap:
        src = infMap[kid]
        idx = user_index[kid]
        if "CASE0" in src:
            status[idx] = 1
        if "PEER" in src:
            status[idx] = 2
            id0 = int(src[5:].split(":")[0])
            idx0 = user_index[id0]
            if status[idx0] == 0:
                status[idx0] = 1
                if print_data_warnings:
                    print("Infecting peer did not have correct status", idx0)
                    
    out = events[events["type"] == "outcome"]
    outMap = pd.Series(out.out.values, index=out.user_id).to_dict()
    for kid in outMap:
        out = outMap[kid]
        idx = user_index[kid]
        if out == "DEAD":
            status[idx] = 3
        if out == "RECOVERED":
            status[idx] = 4
        if out == "VACCINATED":
            status[idx] = 5                    
    
    return status

def get_infection_network(infections, status):
    nvert = len(user_index)
    
    g = Graph(directed=True)
    g.add_vertices(nvert)
    g.add_edges(infections)

    if status:
        g.vs["status"] = status
        g.vs["color"] = [status_color[out] for out in g.vs["status"]]
    
    return g

def plot_network(g, layout, title, fn):
    img_fn = os.path.join(output_folder, fn)
    
    style["layout"] = layout
    p = plot(g, img_fn, **style)
    
    if ".png" in fn and title:
        image = Image.open(img_fn)
        draw = ImageDraw.Draw(image)
        draw.text((10, 760), title, fill='rgb(0, 0, 0)', font=label_font)
        image.save(img_fn)

# https://stackoverflow.com/a/48938464
def hour_rounder(t):
    # Rounds to nearest hour by adding a timedelta hour if minute >= 30
    return (t.replace(second=0, microsecond=0, minute=0, hour=t.hour)
               +timedelta(hours=t.minute//30))

# Load participants and histories

all_users = pd.read_csv(path.join(data_folder, "participants.csv")) 
all_events = pd.read_csv(path.join(data_folder, "histories.csv"))

users = all_users[all_users["sim_id"] == sim_id]

events = all_events[all_events["sim_id"] == sim_id]
events.fillna({'contact_length':0, 'peer_id':-1}, inplace=True)
events["event_start"] = events["time"] - events["contact_length"]/1000
events["event_start"] = events["event_start"].astype(int)

p2pToSim = pd.Series(users.sim_id.values, index=users.p2p_id).to_dict()
p2pToId = pd.Series(users.id.values, index=users.p2p_id).to_dict()
idTop2p = pd.Series(users.p2p_id.values, index=users.id).to_dict()
        
user_index = {}
index_user = {}
idx = 0
for kid in idTop2p:
    user_index[kid] = idx
    index_user[idx] = kid
    idx += 1

# These should return the same value
print(len(users))
print(len(idTop2p))    
print(len(p2pToId))
print(len(user_index))

# Round min and max times to the hour
min_time = min(events['time'])
max_time = max(events['time'])
first_date = hour_rounder(datetime.fromtimestamp(min_time, tz=timezone))
last_date = hour_rounder(datetime.fromtimestamp(max_time, tz=timezone))
min_time = datetime.timestamp(first_date)
max_time = datetime.timestamp(last_date)

print("First event:", first_date)
print("Last event :", last_date)

if time0 and time1:
    print("Start time:", datetime.strptime(time0, '%b %d %Y %I:%M%p'))
    print("End time:", datetime.strptime(time1, '%b %d %Y %I:%M%p'))

print(first_date.tzinfo)

# Infections over time

print("CREATING FRAMES...") 

# How to properly animate an igraph network over time (so nodes change position smoothly from frame to frame):
# http://estebanmoro.org/post/2015-12-21-temporal-networks-with-r-and-igraph-updated/
# https://github.com/emoro/temporal_networks

frame = 0
layout0 = None
tstatus = None

if obs_date0 and obs_date1:
    tmin = datetime.timestamp(obs_date0)
    tmax = datetime.timestamp(obs_date1)
else:
    tmin = min_time
    tmax = max_time

t = tmin
print("FRAME", end =" ")
while t <= tmax:
    t0 = t
    t += time_delta_sec    
    td = datetime.fromtimestamp(t, tz=timezone)    

    condition = events['time'] <= t
    
    tevents = events[condition]
    tstatus = get_node_status(tevents, tstatus)
    infections = get_infection_list(tevents)
    g = get_infection_network(infections, tstatus)
    
    for i in range(0, anim_steps_per_time_delta):
        print(frame, end =" ")
        
        # https://igraph.org/python/api/latest/igraph._igraph.GraphBase.html#layout_fruchterman_reingold        
        layout = g.layout_fruchterman_reingold(niter=fr_niter, start_temp=0.05, grid='nogrid', seed=layout0)
        layout0 = layout.copy()
        
        img_title = td.strftime('%B %d, %I:%M %p')
        img_fn =  "frame-" + str(frame) + "." + frame_format
        plot_network(g, layout, img_title, img_fn)

        frame += 1

print("\nDONE")

print("CREATING THE MOVIE FILE...")

movie_fn = path.join(output_folder, "movie.mp4")
if path.exists(movie_fn):
    os.remove(movie_fn)

cmd_str = "ffmpeg -i " + output_folder + "/frame-%d.png -c:v libx264 -pix_fmt yuv420p " + movie_folder + "/infect-net.mp4"
os.system(cmd_str)

print("DONE")