import os, sys, json
from os import path
from datetime import datetime, timedelta, date
import pytz

from igraph import *

#import openpyxl 

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
status_color = {0: clr.to_hex("cornflowerblue"),  # Susceptible
                1: clr.to_hex("darkorange"),      # Infected (index case)
                2: clr.to_hex("darkorange"),      # Infected (from someone else)
                3: clr.to_hex("darkgrey"),        # Dead 
                4: clr.to_hex("mediumseagreen"),  # Recovered 
                5: clr.to_hex("darkorchid")       # Vaccinated 
               } 

data_folder = path.join(base_folder, "data")
output_folder = path.join(base_folder, "output")
movie_folder = path.join(output_folder, "movies")
if not path.exists(output_folder):
    os.makedirs(output_folder)
output_sir_folder = path.join(base_folder, "output", "charts", "sir")
if not path.exists(output_sir_folder):
    os.makedirs(output_sir_folder)
output_cont_folder = path.join(base_folder, "output", "charts", "contacts")
if not path.exists(output_cont_folder):
    os.makedirs(output_cont_folder)
output_inf_folder = path.join(base_folder, "output", "charts", "infections")
if not path.exists(output_inf_folder):
    os.makedirs(output_inf_folder)       
if not path.exists(movie_folder):
    os.makedirs(movie_folder)

# Print warning messages to the console when parsing data
print_data_warnings = True
    
# Default contact time for transmissions that are missing an associated contact event
def_contact_time = 10
    
frame_format = "png"

# Time delta for plots in seconds
time_delta_sec = 60 * time_step_min

# This hsould match the corresponding parameter inthe infection and contact animations
# so that the animated charts match with them. But for quick renderings, it should be set to 1.
anim_steps_per_time_delta = 30

# Number of ticks in the x axis of epi plots
num_ticks = 10

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

def get_contact_list(events, infections):
    contacts = events[events["type"] == "contact"]

    node0 = contacts.user_id.values
    node1 = contacts.peer_id.values
    length = contacts.contact_length.values

    clist = {}
    for id0, id1, l01 in zip(node0, node1, length):
        n0 = user_index[id0]
        if use_new_id_schema:
            if id1 in user_index:
                n1 = user_index[id1]
            elif print_data_warnings:
                print("Cannot find peer", id1)
        else:
            if id1 in p2pToId:
                n1 = user_index[p2pToId[id1]]
            elif print_data_warnings:
                print("Cannot find peer", id1)
    
        if -1 < n1:
            if n0 < n1:
                p01 = (n0, n1)
            else:
                p01 = (n1, n0)
            if p01 in clist:
                c = clist[p01]
            else: 
                c = 0

            clist[p01] = c + round(l01 / (60 * 1000))
                        
    # Adding contacts from transmissions if they are not registered as contacts already
    for (n0, n1) in infections:
        if n0 < n1:
            p01 = (n0, n1)
        else:
            p01 = (n1, n0)
        if not p01 in clist:
            clist[p01] = def_contact_time
            if print_data_warnings: print("Cannot find contact between", n0, "and", n1)            

    return clist

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
events["event_start"] = events["event_start"].astype(int, errors = 'ignore')
if use_new_id_schema:
    events["peer_id"] = events["peer_id"].astype(int, errors = 'ignore')

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

# Animated charts

frame = 0
tframe = 0
layout0 = None
tstatus = None

if obs_date0 and obs_date1:
    tmin = datetime.timestamp(obs_date0)
    tmax = datetime.timestamp(obs_date1)
    diff_min = (obs_date1 - obs_date0).total_seconds() / 60
else:
    tmin = min_time
    tmax = max_time
    diff_min = (last_date - first_date).total_seconds() / 60

# Calculate label spacing
num_points = diff_min / time_step_min
label_spacing = int(num_points / num_ticks)

nframes = int(((tmax - tmin) / time_delta_sec) * anim_steps_per_time_delta)

series_susceptibles = []
series_infected = []
series_dead = []    
series_recovered = []
series_vaccinated = []
series_contacts = []
series_infections = []
time_index = []
time_ticks = []
tlabels = []

export_susceptibles = []
export_infected = []
export_dead = []    
export_recovered = []
export_vaccinated = []
export_tlabels = []

print("Calculating max number of infections and contacts...", end=" ")
t = tmin
nmaxinf = 0
nmaxcont = 0
while t <= tmax:
    t0 = t    
    t += time_delta_sec
    condition = ((t0 < events['event_start']) & (events['event_start'] <= t)) | ((t0 < events['time']) & (events['time'] <= t))
    tevents = events[condition]
    tinfections = get_infection_list(tevents)
    tcontacts = get_contact_list(tevents, tinfections)
    nmaxinf = max(nmaxinf, len(tinfections))
    nmaxcont = max(nmaxcont, len(tcontacts))    
print("Done")

print("CREATING FRAMES...")
print("FRAME", end =" ")
t = tmin
while t <= tmax:
    t0 = t
    t += time_delta_sec
    td = datetime.fromtimestamp(t, tz=timezone)
        
    # We want to include contact and infection events that either started or ended between t0 and t
    condition = ((t0 < events['event_start']) & (events['event_start'] <= t)) | ((t0 < events['time']) & (events['time'] <= t))
    
    tevents = events[condition]
    tstatus = get_node_status(tevents, tstatus)
    tinfections = get_infection_list(tevents)
    tcontacts = get_contact_list(tevents, tinfections)
    
    nsusceptibles = 0
    ninfected = 0
    ndead = 0    
    nrecovered = 0
    nvaccinated = 0    
    for k in range(0, len(tstatus)):
        if tstatus[k] == 0:
            nsusceptibles += 1
        elif tstatus[k] == 1 or tstatus[k] == 2:
            ninfected += 1  
        elif tstatus[k] == 3:
            ndead += 1
        elif tstatus[k] == 4:
            nrecovered += 1            
        elif tstatus[k] == 5:
            nvaccinated += 1

    ninfections = len(tinfections)            
    ncontacts = len(tcontacts)
            
    ntotal = nsusceptibles + ninfected + ndead + nrecovered + nvaccinated
    # print(nsusceptibles, ninfected, ndead, nrecovered, nvaccinated, ntotal) 

    if tframe % label_spacing == 0:
        tlabels += [td.strftime('%b %d %-I:%M %p')]
        time_ticks += [frame]
    tframe += 1
    
    export_susceptibles.append(nsusceptibles)
    export_infected.append(ninfected)
    export_dead.append(ndead)
    export_recovered.append(nrecovered)
    export_vaccinated.append(nvaccinated)    
    export_tlabels += [td.strftime("%m/%d/%Y %H:%M")]
        
    for i in range(0, anim_steps_per_time_delta):
        print(frame, end =" ")
        series_susceptibles.append(nsusceptibles)
        series_infected.append(ninfected)
        series_dead.append(ndead)
        series_recovered.append(nrecovered)
        series_vaccinated.append(nvaccinated)
        series_contacts.append(ncontacts)
        series_infections.append(ninfections)        
    
        time_index.append(frame)

        # SIR plot
        fig, ax = plt.subplots(figsize=(12,8), facecolor="white")
        plt.ylim([-5, ntotal + 10])
        plt.xlim([-5, nframes + 10])        
        plt.xlabel("Time", labelpad=15, fontsize=15)
        plt.ylabel("Participants", labelpad=15, fontsize=15)
        ax.plot(time_index, series_susceptibles, label="Susceptible", color=status_color[0], lw=2)
        ax.plot(time_index, series_infected, label="Infected", color=status_color[1], lw=2)
        ax.plot(time_index, series_recovered, label="Recovered", color=status_color[4], lw=2)
        ax.plot(time_index, series_vaccinated, label="Vaccinated", color=status_color[5], lw=2)        
        ax.plot(time_index, series_dead, label="Dead", color=status_color[3], lw=2)
        plt.axvline(x=frame, color="dimgray", lw=1)
        plt.xticks(time_ticks, tlabels, rotation=45, horizontalalignment="right")
        plt.legend(loc='upper right')
        plt.tight_layout()
        img_title = td.strftime('%B %d, %I:%M %p')
        img_fn = "frame-" + str(frame) + "." + frame_format
        fig.savefig(os.path.join(output_sir_folder, img_fn))
        plt.close('all')
        
        # Contacts plot
        fig, ax = plt.subplots(figsize=(12,8), facecolor="white")
        plt.ylim([-5, nmaxcont + 10])
        plt.xlim([-5, nframes + 10])        
        plt.xlabel("Time", labelpad=15, fontsize=15)
        plt.ylabel("Number of contacts", labelpad=15, fontsize=15)
        ax.plot(time_index, series_contacts, color="black", lw=2)
        plt.axvline(x=frame, color="dimgray", lw=1)
        plt.xticks(time_ticks, tlabels, rotation=45, horizontalalignment="right")
        plt.tight_layout()
        img_title = td.strftime('%B %d, %I:%M %p')
        img_fn = "frame-" + str(frame) + "." + frame_format
        fig.savefig(os.path.join(output_cont_folder, img_fn))
        plt.close('all')
        
        # Infections plot
        fig, ax = plt.subplots(figsize=(12,8), facecolor="white")
        plt.ylim([-5, nmaxinf + 10])
        plt.xlim([-5, nframes + 10])        
        plt.xlabel("Time", labelpad=15, fontsize=15)
        plt.ylabel("Number of infections", labelpad=15, fontsize=15)
        ax.plot(time_index, series_infections, color=status_color[1], lw=2)
        plt.axvline(x=frame, color="dimgray", lw=1)
        plt.xticks(time_ticks, tlabels, rotation=45, horizontalalignment="right")
        plt.tight_layout()
        img_title = td.strftime('%B %d, %I:%M %p')
        img_fn = "frame-" + str(frame) + "." + frame_format
        fig.savefig(os.path.join(output_inf_folder, img_fn))
        plt.close('all') 

        frame += 1

print("\nDONE")

# Saving data file
df = pd.DataFrame({"Time": export_tlabels,                    
                   "Susceptible": export_susceptibles, "Infected": export_infected, "Dead": export_dead, "Recovered": export_recovered, "Vaccinated": export_vaccinated})
df.to_excel(os.path.join(output_folder, "epi-data.xlsx"), index=False)

print("CREATING THE MOVIE FILES...")

def make_movie(in_folder, out_folder, fn):
    movie_fn = path.join(out_folder, fn)
    if path.exists(movie_fn):
        os.remove(movie_fn)
    cmd_str = "ffmpeg -i " + in_folder + "/frame-%d.png -c:v libx264 -pix_fmt yuv420p " + out_folder + "/" + fn
    os.system(cmd_str)

make_movie(output_sir_folder, movie_folder, "counts-sir.mp4")
make_movie(output_cont_folder, movie_folder, "counts-cont.mp4")
make_movie(output_inf_folder, movie_folder, "counts-inf.mp4")

print("DONE")

# R effective over time

# Might need a scaling larger than 1 to capture more events for an accurate estimation of Reff
scale = 5

delta = time_delta_sec * scale
spacing = int(label_spacing / scale)

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
r_mean_values = []
r_std_values = []
tlabels = []
time_ticks = []
index = []
while t <= tmax:
    t0 = t
    t += delta
    td = datetime.fromtimestamp(t, tz=timezone)    

    # We want to include contact and infection events that either started or ended between t0 and t
    condition = ((t0 < events['event_start']) & (events['event_start'] <= t)) | ((t0 < events['time']) & (events['time'] <= t))
    
    tevents = events[condition]
    tstatus = get_node_status(tevents, tstatus)
    tinfections = get_infection_list(tevents)
    g = get_infection_network(tinfections, tstatus)
    
    # Getting all nodes with at least one edge
    r_values = []    
    sel = g.vs(_degree_gt=0)

    for v in sel:
        nout = v.degree(mode=OUT)
        r_values += [nout]

    if r_values:
        r_mean = np.mean(r_values)
        r_std = np.std(r_values)
        r_mean_values += [r_mean]
        r_std_values += [r_std]
    else:
        r_mean_values += [0]
        r_std_values += [0]
        
    if frame % spacing == 0: 
        tlabels += [td.strftime('%b %d %-I:%M %p')]
        time_ticks += [frame]
    time_index += [frame]
    frame += 1

mu = np.array(r_mean_values)
sigma = np.array(r_std_values)
time = np.arange(len(sigma))
    
fig, ax = plt.subplots(figsize=(12,8))
plt.xlabel("Time", labelpad=15, fontsize=15)
plt.ylabel("R", labelpad=15, fontsize=15)
ax.plot(time, mu, lw=2, label='R effective', color='darkorange')
ax.fill_between(time, mu+sigma, mu-sigma, facecolor='darkorange', alpha=0.5)

plt.xticks(time_ticks, tlabels, rotation=45, horizontalalignment="right")

plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
fig.savefig(os.path.join(output_folder, "r-effective.pdf"))