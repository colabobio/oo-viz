# Operation Outbreak visualization scripts

This repo contains Jupyter notebooks and other tools to visualize OO data:

* plot_contacts.ipynb: Generates an animated force-directed network graph showing the contacts and infections during the simulations
* plot_infections.ipynb: Generates an animated force-directed network graph showing the infection chains as the appear and grow during the simulation
* plot_charts.ipynb: Generates animated 2D charts showing the number of susceptible, infected and removed players over time, number of contacts over time, and number of new cases over time.

In all cases, the notebooks generate a movie file that can be used to play the animations outside the notebook environment.

## Dependencies

The notebook uses some Python librariews for plotting:

* [python-igraph](https://igraph.org/python/)
* [seaborn](https://seaborn.pydata.org/)

You can install the requirements by running 

```pip install -r requirements.txt```

## Resources to visualize epidemiological data

The [Outbreakteach R project](https://github.com/mrc-ide/outbreakteachR) is a good reference for epi visualization code.

Some options for graph/network visualization:

* [igraph](https://igraph.org)
* [networkx](https://networkx.github.io/)
* [graph-tool](https://graph-tool.skewed.de/)
* [Performance comparison](https://www.timlrx.com/2019/05/05/benchmark-of-popular-graph-network-packages/) between those and a few other packages.
