# 09: Network plot

**What to build:** `plot_network(de)` draws the evidence network — drugs as nodes, comparative estimates as edges — built on igraph/ggraph, with optional visual encodings: edge thickness = precision, node shape = presence of an absolute anchor, node size = number of comparisons, drug-class grouping, edge labels showing observed effects.

**Blocked by:** 02 (Direct-effect network object)

**Status:** ready-for-agent

- [ ] Works on a bare network object — no fit required
- [ ] Optional encodings can be toggled, with sensible defaults
- [ ] Graph logic lives with the network object / igraph, not inside plotting code
- [ ] Renders correctly for multi-component networks
