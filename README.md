# O2 visualization tools

This repo contains some scripts and other tools to visualize O2 data.

To visualize the log data from a simulation, run the plot_data Jupyter notebook

## Dependencies

The notebook uses some Python librariews for plotting:

* [python-igraph](https://igraph.org/python/)
* [seaborn](https://seaborn.pydata.org/)

You can install the requirements by running 

```pip install -r requirements.txt```

## Resources to visualize epidemiological data

The [Outbreakteach R project](https://github.com/mrc-ide/outbreakteachR) is a good reference for epi visualization code.

Some options for graph/network visualization:

[igraph](https://igraph.org)

[networkx](https://networkx.github.io/)

[graph-tool](https://graph-tool.skewed.de/)

[Performance comparison](https://www.timlrx.com/2019/05/05/benchmark-of-popular-graph-network-packages/) between those and a few other packages.