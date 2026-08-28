# 09: Network plot

**What to build:** `plot_network(de)` draws the evidence network — drugs as nodes, comparative estimates as edges — built on igraph/ggraph, with optional visual encodings: edge thickness = precision, node shape = presence of an absolute anchor, node size = number of comparisons, drug-class grouping, edge labels showing observed effects.

**Blocked by:** 02 (Direct-effect network object)

**Status:** done

- [x] Works on a bare network object — no fit required
- [x] Optional encodings can be toggled, with sensible defaults
- [x] Graph logic lives with the network object / igraph, not inside plotting code
- [x] Renders correctly for multi-component networks

## Comments

Done in R/plot-network.R on ggraph (Suggests) with the deterministic
"stress" layout, which also packs disconnected components.
network_node_data() (R/network.R) derives per-drug comparison counts
and anchor flags — graph logic stays with the network module; the
constructor's igraph graph already carries study/estimate/std_error
edge attributes. Encodings, each toggleable: edge width = precision
(default on), node size = comparison count (on), node shape = anchor
presence (triangle/circle, on), edge labels = observed effects (off),
node colour = drug class via an explicit `drug_classes` data frame
(drug, class) — per-comparison drug_class metadata is ambiguous about
which drug it describes, so classes are supplied per drug. Repeat
comparisons draw as fanned parallel edges. Tests force full rendering
through ggplot_build(), including multi-component and repeat-edge
networks. Local apt ggraph 2.1.0 emits ggplot2 size-deprecation
warnings from its own internals; current CRAN ggraph in CI does not.
