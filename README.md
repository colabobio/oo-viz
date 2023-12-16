# Operation Outbreak visualization scripts

This repo contains Jupyter notebooks and other tools to visualize OO data:

* plot_contacts.ipynb: Generates an animated force-directed network graph showing the contacts and infections during the simulations
* plot_infections.ipynb: Generates an animated force-directed network graph showing the infection chains as the appear and grow during the simulation
* plot_charts.ipynb: Generates animated 2D charts showing the number of susceptible, infected and removed players over time, number of contacts over time, and number of new cases over time.
* plot_sequences.ipynb: Generates the sequences associates to the transmissions by applying the provided mutations to the reference sequence, and plots the network of sequence transmissions and the phylogenetic tree based on the alignment of the sequences.

In all cases, the notebooks generate a movie file that can be used to play the animations outside the notebook environment.

To visualiza the data for very large simulationd, which may hang the notebook, use the corresponding python scripts that can be run from the command line using the simulaton properties in json format:

* python plot_contacts.py simulations/ootest/sim.json
* python plot_infections.py simulations/ootest/sim.json
* python plot_charts.py simulations/ootest/sim.json
* python plot_sequences.py simulations/ootest/sim.json

The json file shoud have the following format:

```
{
    "title": "SIM TITLE",
    "base_folder": "./simulations/ootest",
    "sim_id": 111,
    "pathogen_id": 11,
    "sim_tz": "America/New York",
    "time0": "Nov 20 2023 9:00AM",
    "time1": "Dec 4 2023 12:00PM",
    "time_step_min": 30,
    "use_new_id_schema": true
}
```

## Dependencies

The notebook uses some Python librariews for plotting:

* [python-igraph](https://igraph.org/python/)
* [seaborn](https://seaborn.pydata.org/)

When using conda, it is recommended to create a separate environment to load all the dependencies.

Then, install cairo and ffmpeg with brew:

```
brew install cairo
brew install ffmpeg
```

Then all the Python dependencies can be installed by running:

```pip install -r requirements.txt```

If ffmpeg is crashing with an error like this:

```
$ ffmpeg -i simulations/wku23/output/charts/sir/frame-%d.png -c:v libx264 -pix_fmt yuv420p simulations/wku23/output/movies/counts-sir.mp4
dyld[90963]: Library not loaded: /opt/homebrew/opt/libunistring/lib/libunistring.2.dylib
  Referenced from: <DC03D607-E12F-36EF-B99B-F9FAE69BDA4E> /opt/homebrew/Cellar/libidn2/2.3.4/lib/libidn2.0.dylib
  Reason: tried: '/opt/homebrew/opt/libunistring/lib/libunistring.2.dylib' (no such file), '/System/Volumes/Preboot/Cryptexes/OS/opt/homebrew/opt/libunistring/lib/libunistring.2.dylib' (no such file), '/opt/homebrew/opt/libunistring/lib/libunistring.2.dylib' (no such file), '/usr/local/lib/libunistring.2.dylib' (no such file), '/usr/lib/libunistring.2.dylib' (no such file, not in dyld cache), '/opt/homebrew/Cellar/libunistring/1.1/lib/libunistring.2.dylib' (no such file), '/System/Volumes/Preboot/Cryptexes/OS/opt/homebrew/Cellar/libunistring/1.1/lib/libunistring.2.dylib' (no such file), '/opt/homebrew/Cellar/libunistring/1.1/lib/libunistring.2.dylib' (no such file), '/usr/local/lib/libunistring.2.dylib' (no such file), '/usr/l
```

create the mising dependency by linking to the correct unistring library file, e.g.:

```ln -s /opt/homebrew/opt/libunistring/lib/libunistring.5.dylib  /opt/homebrew/opt/libunistring/lib/libunistring.2.dylib```

## Resources to visualize epidemiological data

The [Outbreakteach R project](https://github.com/mrc-ide/outbreakteachR) is a good reference for epi visualization code.

Some options for graph/network visualization:

* [igraph](https://igraph.org)
* [networkx](https://networkx.github.io/)
* [graph-tool](https://graph-tool.skewed.de/)
* [Performance comparison](https://www.timlrx.com/2019/05/05/benchmark-of-popular-graph-network-packages/) between those and a few other packages.
