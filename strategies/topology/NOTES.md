# Topology strategy: notes

All numbers below were computed from the delivered files (`topology_<algorithm>/out/<model>/`) and from `strategy.py` at seed 1 unless a seed range is stated. All six panel models were run with all three algorithms; nothing was dropped or subsampled. Sizes count nodes (species and reactions), currency species included in their home group.

Layout: `strategy.py`, `RULES.md` and this file are shared; each variant folder `topology_<algorithm>/` holds `out/<model>/` (`groups.tsv`, `groups_meta.tsv`, `<model>_grouped.xml`) and `stats/<model>.json` (the per-run numbers behind the tables, written by `--stats`).

**Recommended variant: `topology_leiden`.** Reasons are in section 3.

## 1. Per model

### brenda_merged

396 nodes. Degree bound of R2: 13. Candidates: 15; cut as currency species by R3: 15. The currency step is the same for the three variants.

| candidate | degree | share of its reactions stranded | result |
|---|---|---|---|
| H2O (`M_mt_1b1223fff4`) | 65 | 0.00 | currency (cut, cloned) |
| ATP (`M_mt_82e3f482bf`) | 43 | 0.00 | currency (cut, cloned) |
| ADP (`M_mt_7133524672`) | 42 | 0.02 | currency (cut, cloned) |
| phosphate (`M_mt_3577006728`) | 36 | 0.03 | currency (cut, cloned) |
| NH3 (`M_mt_a582428d9f`) | 24 | 0.00 | currency (cut, cloned) |
| NADH (`M_mt_19acfa327b`) | 21 | 0.00 | currency (cut, cloned) |
| NAD+ (`M_mt_62968a0d8b`) | 21 | 0.00 | currency (cut, cloned) |
| H+ (`M_mt_d78fb7f9bf`) | 20 | 0.00 | currency (cut, cloned) |
| CO2 (`M_mt_9cf23c7e50`) | 19 | 0.00 | currency (cut, cloned) |
| L-glutamate (`M_mt_219b3a1a3f`) | 17 | 0.06 | currency (cut, cloned) |
| NADP+ (`M_mt_410bd09215`) | 16 | 0.00 | currency (cut, cloned) |
| NADPH (`M_mt_942f98c2ec`) | 16 | 0.00 | currency (cut, cloned) |
| AMP (`M_mt_7fdd97c6e1`) | 14 | 0.07 | currency (cut, cloned) |
| diphosphate (`M_mt_bdcc7261de`) | 13 | 0.00 | currency (cut, cloned) |
| 2-oxoglutarate (`M_mt_ea23b49b7a`) | 13 | 0.15 | currency (cut, cloned) |

| | leiden | infomap | walktrap |
|---|---|---|---|
| groups L1 / L2 / L3 | 12 / 0 / 0 | 12 / 8 / 0 | 12 / 4 / 0 |
| L1 size min / median / max | 18 / 29 / 56 | 17 / 29 / 64 | 12 / 28 / 62 |
| L2 size min / median / max | - | 5 / 8 / 13 | 10 / 15 / 22 |
| L3 size min / median / max | - | - | - |
| leaf groups (no children): number | 12 | 19 | 15 |
| leaf size min / median / max | 18 / 29 / 56 | 5 / 21 / 48 | 10 / 22 / 61 |
| largest L1 group, share of all nodes | 14.1% | 16.2% | 15.7% |
| leaf groups below 5 nodes / above 60 nodes | 0 / 0 | 0 / 0 | 0 / 1 |
| unnamed groups (share) | 0 of 12 (0%) | 0 of 20 (0%) | 0 of 16 (0%) |
| name sources | member 12 | ec 1, member 19 | ec 2, member 14 |
| duplicate labels | 0 | 0 | 0 |
| modules returned by the algorithm (R4, before R5-R6) | 22 | 47 | 31 |
| modularity: R4 output / L1 / leaf groups | 0.795 / 0.789 / 0.789 | 0.741 / 0.783 / 0.761 | 0.756 / 0.761 / 0.763 |
| nodes placed by the algorithm itself (rule R4 or R7) | 83% | 48% | 71% |
| nodes per rule | R4 328, R5b 8, R6 45, R8 15 | R4 131, R5a 11, R5b 8, R6 174, R7 57, R8 15 | R4 231, R5a 20, R5b 8, R6 71, R7 51, R8 15 |
| currency species with clones / clones in total | 15 / 107 | 15 / 134 | 15 / 108 |
| leaf groups lying in one compartment (`loc_space`) | 100% | 100% | 100% |
| nodes in the majority compartment of their leaf group | 100% | 100% | 100% |

Agreement between the variants (adjusted Rand index / normalised mutual information over all nodes):

| pair | L1: ARI / NMI | leaf groups: ARI / NMI |
|---|---|---|
| leiden vs infomap | 0.76 / 0.85 | 0.69 / 0.83 |
| leiden vs walktrap | 0.59 / 0.73 | 0.61 / 0.75 |
| infomap vs walktrap | 0.51 / 0.70 | 0.41 / 0.69 |

Stability (seeds 1 to 10, adjusted Rand index over the 45 pairs of runs; the delivered files are seed 1):

| | L1: mean (min) | leaf groups: mean (min) | modules from R4, range | modularity of R4 output, range |
|---|---|---|---|---|
| leiden | 0.64 (0.47) | 0.62 (0.47) | 20 to 24 | 0.789 to 0.795 |
| infomap | 0.69 (0.52) | 0.66 (0.46) | 46 to 52 | 0.725 to 0.743 |
| walktrap | no random step: identical for every seed (checked, seeds 1 and 7) | same | - | - |

Walktrap instead depends on the walk length. Against the delivered 4 steps: 3 steps ARI 0.68 (L1) and 0.74 (leaf groups); 5 steps ARI 0.65 (L1) and 0.68 (leaf groups).

Three example groups of `topology_leiden` (the largest, the median and the smallest leaf group):

- **M01** (largest), label "IMP" (`name_source` member, `name_xref` `M_mt_15d5819dc2`), 56 nodes: 22 species, 34 reactions. Species: 5-phospho-alpha-D-ribose 1-diphosphate, ADP-ribose, AMP, D-ribose, D-ribose 5-phosphate, GMP, H2O, IMP, Inosine, NH3, and others. Reactions: BR149; BR16042; BR16044; BR19230, and others.
- **M08** (median), label "formaldehyde" (`name_source` member, `name_xref` `M_mt_224abb736c`), 29 nodes: 16 species, 13 reactions. Species: 6,8-dihydroxypurine, 7-methylxanthine, 8-hydroxypurine, Caffeine, D-arabino-hex-3-ulose 6-phosphate, H+, NAD+, NADH, O2, Theobromine, and others. Reactions: BR48027; BR48972; BR48975; BR50577-a, and others.
- **M12** (smallest), label "DL-glyceraldehyde" (`name_source` member, `name_xref` `M_mt_7476acb4c8`), 18 nodes: 9 species, 9 reactions. Species: (R)-2-hydroxyglutaryl-1-CoA, D-glycerate, DL-glyceraldehyde, acceptor, butyryl-CoA, crotonyl-CoA, glutaconyl-CoA, oxidized acceptor, reduced acceptor. Reactions: BR48991; BR50094; BR51549; BS370908, and others.

### e_coli_core

167 nodes. Degree bound of R2: 11. Candidates: 9; cut as currency species by R3: 9. The currency step is the same for the three variants.

| candidate | degree | share of its reactions stranded | result |
|---|---|---|---|
| H+ (`M_h_c`) | 35 | 0.00 | currency (cut, cloned) |
| H2O H2O (`M_h2o_c`) | 18 | 0.06 | currency (cut, cloned) |
| H+ (`M_h_e`) | 17 | 0.06 | currency (cut, cloned) |
| ATP C10H12N5O13P3 (`M_atp_c`) | 13 | 0.00 | currency (cut, cloned) |
| ADP C10H12N5O10P2 (`M_adp_c`) | 12 | 0.00 | currency (cut, cloned) |
| Nicotinamide adenine dinucleotide (`M_nad_c`) | 12 | 0.00 | currency (cut, cloned) |
| Nicotinamide adenine dinucleotide - reduced (`M_nadh_c`) | 12 | 0.00 | currency (cut, cloned) |
| Phosphate (`M_pi_c`) | 12 | 0.25 | currency (cut, cloned) |
| Pyruvate (`M_pyr_c`) | 11 | 0.18 | currency (cut, cloned) |

| | leiden | infomap | walktrap |
|---|---|---|---|
| groups L1 / L2 / L3 | 11 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| L1 size min / median / max | 5 / 9 / 34 | 5 / 12 / 37 | 5 / 8 / 53 |
| L2 size min / median / max | - | - | - |
| L3 size min / median / max | - | - | - |
| leaf groups (no children): number | 11 | 12 | 12 |
| leaf size min / median / max | 5 / 9 / 34 | 5 / 12 / 37 | 5 / 8 / 53 |
| largest L1 group, share of all nodes | 20.4% | 22.2% | 31.7% |
| leaf groups below 5 nodes / above 60 nodes | 0 / 0 | 0 / 0 | 0 / 0 |
| unnamed groups (share) | 0 of 11 (0%) | 0 of 12 (0%) | 0 of 12 (0%) |
| name sources | member 11 | ec 1, member 11 | ec 1, member 11 |
| duplicate labels | 0 | 0 | 0 |
| modules returned by the algorithm (R4, before R5-R6) | 20 | 31 | 37 |
| modularity: R4 output / L1 / leaf groups | 0.703 / 0.692 / 0.692 | 0.668 / 0.688 / 0.688 | 0.627 / 0.656 / 0.656 |
| nodes placed by the algorithm itself (rule R4 or R7) | 82% | 56% | 54% |
| nodes per rule | R4 137, R5a 9, R5b 12, R8 9 | R4 93, R5a 31, R5b 12, R6 22, R8 9 | R4 90, R5a 50, R5b 12, R6 6, R8 9 |
| currency species with clones / clones in total | 9 / 52 | 9 / 57 | 9 / 50 |
| leaf groups lying in one compartment (`loc_space`) | 18% | 25% | 25% |
| nodes in the majority compartment of their leaf group | 76% | 76% | 76% |

Agreement between the variants (adjusted Rand index / normalised mutual information over all nodes):

| pair | L1: ARI / NMI | leaf groups: ARI / NMI |
|---|---|---|
| leiden vs infomap | 0.62 / 0.81 | 0.62 / 0.81 |
| leiden vs walktrap | 0.48 / 0.73 | 0.48 / 0.73 |
| infomap vs walktrap | 0.62 / 0.83 | 0.62 / 0.83 |

Stability (seeds 1 to 10, adjusted Rand index over the 45 pairs of runs; the delivered files are seed 1):

| | L1: mean (min) | leaf groups: mean (min) | modules from R4, range | modularity of R4 output, range |
|---|---|---|---|---|
| leiden | 0.73 (0.46) | 0.73 (0.46) | 17 to 20 | 0.702 to 0.708 |
| infomap | 0.60 (0.40) | 0.60 (0.40) | 30 to 32 | 0.654 to 0.679 |
| walktrap | no random step: identical for every seed (checked, seeds 1 and 7) | same | - | - |

Walktrap instead depends on the walk length. Against the delivered 4 steps: 3 steps ARI 0.93 (L1) and 0.93 (leaf groups); 5 steps ARI 0.93 (L1) and 0.93 (leaf groups).

Three example groups of `topology_leiden` (the largest, the median and the smallest leaf group):

- **M01** (largest), label "Nicotinamide adenine dinucleotide phosphate" (`name_source` member, `name_xref` `M_nadp_c`), 34 nodes: 15 species, 19 reactions. Species: 2-Oxoglutarate, Ammonium, H+, L-Glutamate, L-Glutamine, Nicotinamide adenine dinucleotide, Nicotinamide adenine dinucleotide - reduced, Nicotinamide adenine dinucleotide phosphate, Nicotinamide adenine dinucleotide phosphate - reduced, Phosphate. Reactions: 2 oxoglutarate reversible transport via symport; 2-Oxoglutarate exchange; Ammonia exchange; Ammonia reversible transport, and others.
- **M06** (median), label "D-Glucose 6-phosphate" (`name_source` member, `name_xref` `M_g6p_c`), 9 nodes: 4 species, 5 reactions. Species: 6-Phospho-D-gluconate, 6-phospho-D-glucono-1,5-lactone, D-Glucose, D-Glucose 6-phosphate. Reactions: 6-phosphogluconolactonase; D-Glucose exchange; D-glucose transport via PEP:Pyr PTS; Glucose 6-phosphate dehydrogenase, and others.
- **M11** (smallest), label "D-Lactate" (`name_source` member, `name_xref` `M_lac__D_c`), 5 nodes: 2 species, 3 reactions. Species: D-Lactate. Reactions: D lactate transport via proton symport; D-lactate dehydrogenase; D-lactate exchange.

### R-HSA-77289

148 nodes. Degree bound of R2: 8. Candidates: 9; cut as currency species by R3: 9. The currency step is the same for the three variants.

| candidate | degree | share of its reactions stranded | result |
|---|---|---|---|
| H+ [mitochondrial matrix] (`species_113529`) | 20 | 0.00 | currency (cut, cloned) |
| CoA-SH [mitochondrial matrix] (`species_29374`) | 13 | 0.08 | currency (cut, cloned) |
| FAD [mitochondrial matrix] (`species_113596`) | 12 | 0.00 | currency (cut, cloned) |
| FADH2 [mitochondrial matrix] (`species_31649`) | 12 | 0.17 | currency (cut, cloned) |
| Trifunctional Protein [mitochondrial matrix] (`species_77267`) | 12 | 0.00 | currency (cut, cloned) |
| H2O [mitochondrial matrix] (`species_113521`) | 10 | 0.10 | currency (cut, cloned) |
| NAD+ [mitochondrial matrix] (`species_113526`) | 10 | 0.00 | currency (cut, cloned) |
| NADH [mitochondrial matrix] (`species_29362`) | 10 | 0.00 | currency (cut, cloned) |
| Ac-CoA [mitochondrial matrix] (`species_113559`) | 9 | 0.11 | currency (cut, cloned) |

| | leiden | infomap | walktrap |
|---|---|---|---|
| groups L1 / L2 / L3 | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| L1 size min / median / max | 4 / 10 / 34 | 4 / 12 / 24 | 4 / 11 / 28 |
| L2 size min / median / max | - | - | - |
| L3 size min / median / max | - | - | - |
| leaf groups (no children): number | 12 | 12 | 12 |
| leaf size min / median / max | 4 / 10 / 34 | 4 / 12 / 24 | 4 / 11 / 28 |
| largest L1 group, share of all nodes | 23.0% | 16.2% | 18.9% |
| leaf groups below 5 nodes / above 60 nodes | 2 / 0 | 2 / 0 | 2 / 0 |
| unnamed groups (share) | 0 of 12 (0%) | 0 of 12 (0%) | 0 of 12 (0%) |
| name sources | ec 1, go 1, member 10 | ec 2, member 10 | ec 1, member 11 |
| duplicate labels | 0 | 0 | 0 |
| modules returned by the algorithm (R4, before R5-R6) | 16 | 24 | 22 |
| modularity: R4 output / L1 / leaf groups | 0.799 / 0.775 / 0.775 | 0.756 / 0.762 / 0.762 | 0.762 / 0.767 / 0.767 |
| nodes placed by the algorithm itself (rule R4 or R7) | 78% | 48% | 59% |
| nodes per rule | R4 115, R5b 16, R5c 8, R8 9 | R4 71, R5a 4, R5b 16, R5c 8, R6 40, R8 9 | R4 87, R5a 11, R5b 16, R5c 8, R6 17, R8 9 |
| currency species with clones / clones in total | 9 / 47 | 9 / 37 | 9 / 44 |
| leaf groups lying in one compartment (`loc_space`) | 100% | 100% | 100% |
| nodes in the majority compartment of their leaf group | 100% | 100% | 100% |

Agreement between the variants (adjusted Rand index / normalised mutual information over all nodes):

| pair | L1: ARI / NMI | leaf groups: ARI / NMI |
|---|---|---|
| leiden vs infomap | 0.43 / 0.75 | 0.43 / 0.75 |
| leiden vs walktrap | 0.58 / 0.81 | 0.58 / 0.81 |
| infomap vs walktrap | 0.53 / 0.77 | 0.53 / 0.77 |

Stability (seeds 1 to 10, adjusted Rand index over the 45 pairs of runs; the delivered files are seed 1):

| | L1: mean (min) | leaf groups: mean (min) | modules from R4, range | modularity of R4 output, range |
|---|---|---|---|---|
| leiden | 0.64 (0.42) | 0.64 (0.42) | 15 to 16 | 0.799 to 0.804 |
| infomap | 0.64 (0.49) | 0.64 (0.49) | 22 to 25 | 0.751 to 0.775 |
| walktrap | no random step: identical for every seed (checked, seeds 1 and 7) | same | - | - |

Walktrap instead depends on the walk length. Against the delivered 4 steps: 3 steps ARI 0.48 (L1) and 0.48 (leaf groups); 5 steps ARI 0.73 (L1) and 0.73 (leaf groups).

Three example groups of `topology_leiden` (the largest, the median and the smallest leaf group):

- **M01** (largest), label "HADH dimer [mitochondrial matrix]" (`name_source` member, `name_xref` `species_71052`), 34 nodes: 23 species, 11 reactions. Species: (S)-Hydroxyhexanoyl-CoA [mitochondrial matrix], 3-Oxohexanoyl-CoA [mitochondrial matrix], 3-Oxooctanoyl-CoA [mitochondrial matrix], 3HB-CoA [mitochondrial matrix], 3OA-CoA [mitochondrial matrix], ACA-CoA [mitochondrial matrix], ACAA2 tetramer [mitochondrial matrix], ACADS tetramer [mitochondrial matrix], ACOT2,9, THEM4,5 dimer [mitochondrial matrix], Ac-CoA [mitochondrial matrix], and others. Reactions: (S)-Hydroxybutanoyl-CoA+NAD => Acetoacetyl-CoA+NADH+H; (S)-Hydroxyhexanoyl-CoA+NAD => 3-Oxohexanoyl-CoA+NADH+H; (S)-Hydroxyoctanoyl-CoA+NAD => 3-Oxooctanoyl-CoA+NADH+H; 3-Oxohexanoyl-CoA+CoA-SH => Butanoyl-CoA + acetyl-CoA, and others.
- **M07** (median), label "EC 6.2.1" (`name_source` ec, `name_xref` `ec-code:6.2.1`), 10 nodes: 8 species, 2 reactions. Species: ACSF2 [mitochondrial matrix], ACSM3,ACSM6 [mitochondrial matrix], AMP [mitochondrial matrix], ATP [mitochondrial matrix], BUT [mitochondrial matrix], MCFA [mitochondrial matrix], MCFA-CoA [mitochondrial matrix], PPi [mitochondrial matrix]. Reactions: ACSF2 ligates CoA-SH to MCFA; ACSM3,ACSM6 ligate CoA to BUT.
- **M12** (smallest), label "ACBD6 binds LCFA-CoA" (`name_source` member, `name_xref` `reaction_8848250`), 4 nodes: 3 species, 1 reactions. Species: ACBD6 [cytosol], ACBD6:LCFA-CoA [cytosol], LCFA-CoA [cytosol]. Reactions: ACBD6 binds LCFA-CoA.

### skm_pss

1817 nodes. Degree bound of R2: 9. Candidates: 19; cut as currency species by R3: 13. The currency step is the same for the three variants.

| candidate | degree | share of its reactions stranded | result |
|---|---|---|---|
| VPg (`s_VPg_cyt_p`) | 53 | 0.83 | stays in the graph (hub) |
| SARD1 (`s_SARD1_nuc_pa`) | 31 | 0.00 | currency (cut, cloned) |
| CBP60G (`s_CBP60G_nuc_pa`) | 30 | 0.53 | stays in the graph (hub) |
| CDPK (`s_CDPK_cyt_pa`) | 25 | 0.12 | currency (cut, cloned) |
| MYC2 (`s_MYC2_nuc_pa`) | 21 | 0.24 | currency (cut, cloned) |
| Ca2+ (`s_Ca2_cyt_m`) | 20 | 0.15 | currency (cut, cloned) |
| HC-Pro (`s_HCPro_cyt_pa`) | 18 | 0.44 | currency (cut, cloned) |
| SNRK2 (`s_SNRK2_cyt_pa`) | 18 | 0.06 | currency (cut, cloned) |
| SA (`s_SA_cyt_m`) | 16 | 0.19 | currency (cut, cloned) |
| ORA59 (`s_ORA59_nuc_pa`) | 15 | 0.67 | stays in the graph (hub) |
| MPK3,6 (`s_MPK36_cyt_pa`) | 14 | 0.21 | currency (cut, cloned) |
| ICS (`s_ICS_chl_pa`) | 13 | 0.15 | currency (cut, cloned) |
| ROS (`s_ROS_cyt_m`) | 12 | 0.17 | currency (cut, cloned) |
| SNRK2 (`s_SNRK2_cyt_p`) | 12 | 0.83 | stays in the graph (hub) |
| IAA (`s_IAA_cyt_m`) | 10 | 0.80 | stays in the graph (hub) |
| NPR1 (`s_NPR1_nuc_pa`) | 10 | 0.50 | stays in the graph (hub) |
| BIK1 (`s_BIK1_cyt_pa`) | 9 | 0.33 | currency (cut, cloned) |
| JAZ (`s_JAZ_nuc_pa`) | 9 | 0.11 | currency (cut, cloned) |
| potyvirus (`s_potyvirus_cyt_f`) | 9 | 0.22 | currency (cut, cloned) |

| | leiden | infomap | walktrap |
|---|---|---|---|
| groups L1 / L2 / L3 | 12 / 96 / 0 | 12 / 96 / 0 | 12 / 95 / 6 |
| L1 size min / median / max | 4 / 197 / 295 | 4 / 180 / 298 | 4 / 193 / 359 |
| L2 size min / median / max | 5 / 16 / 133 | 5 / 18 / 133 | 5 / 14 / 155 |
| L3 size min / median / max | - | - | 7 / 13 / 146 |
| leaf groups (no children): number | 100 | 100 | 103 |
| leaf size min / median / max | 4 / 16 / 133 | 4 / 16 / 133 | 4 / 13 / 146 |
| largest L1 group, share of all nodes | 16.2% | 16.4% | 19.8% |
| leaf groups below 5 nodes / above 60 nodes | 1 / 1 | 1 / 1 | 1 / 1 |
| unnamed groups (share) | 0 of 108 (0%) | 0 of 108 (0%) | 0 of 113 (0%) |
| name sources | member 96, notes 12 | ec 1, member 96, notes 11 | member 102, notes 11 |
| duplicate labels | 0 | 0 | 0 |
| modules returned by the algorithm (R4, before R5-R6) | 68 | 261 | 133 |
| modularity: R4 output / L1 / leaf groups | 0.925 / 0.850 / 0.894 | 0.834 / 0.847 / 0.889 | 0.879 / 0.835 / 0.883 |
| nodes placed by the algorithm itself (rule R4 or R7) | 83% | 45% | 66% |
| nodes per rule | R4 24, R5a 92, R5b 55, R5c 4, R6 137, R7 1492, R8 13 | R4 24, R5a 257, R5b 58, R5c 4, R6 675, R7 786, R8 13 | R4 20, R5a 309, R5b 58, R5c 4, R6 224, R7 1189, R8 13 |
| currency species with clones / clones in total | 13 / 96 | 13 / 92 | 13 / 90 |
| leaf groups lying in one compartment (`loc_space`) | 36% | 42% | 40% |
| nodes in the majority compartment of their leaf group | 83% | 83% | 81% |

Agreement between the variants (adjusted Rand index / normalised mutual information over all nodes):

| pair | L1: ARI / NMI | leaf groups: ARI / NMI |
|---|---|---|
| leiden vs infomap | 0.34 / 0.51 | 0.80 / 0.91 |
| leiden vs walktrap | 0.26 / 0.45 | 0.71 / 0.89 |
| infomap vs walktrap | 0.23 / 0.42 | 0.69 / 0.88 |

Stability (seeds 1 to 10, adjusted Rand index over the 45 pairs of runs; the delivered files are seed 1):

| | L1: mean (min) | leaf groups: mean (min) | modules from R4, range | modularity of R4 output, range |
|---|---|---|---|---|
| leiden | 0.42 (0.29) | 0.83 (0.81) | 67 to 71 | 0.925 to 0.926 |
| infomap | 0.42 (0.28) | 0.82 (0.78) | 255 to 261 | 0.833 to 0.837 |
| walktrap | no random step: identical for every seed (checked, seeds 1 and 7) | same | - | - |

Walktrap instead depends on the walk length. Against the delivered 4 steps: 3 steps ARI 0.28 (L1) and 0.74 (leaf groups); 5 steps ARI 0.34 (L1) and 0.75 (leaf groups).

Three example groups of `topology_leiden` (the largest, the median and the smallest leaf group):

- **M01.01** (largest), label "6K2|VPg" (`name_source` member, `name_xref` `s_6K2VPg_cyt_c`), 133 nodes: 88 species, 45 reactions. Species: 6K2, 6K2|VPg, ABCF1, ABCF1|VPg, AT3G17020, AT3G17020|VPg, ATPC1, ATPC1|VPg, BAS1, BAS1|VPg, and others. The reactions of this model have no names.
- **M04.08** (median), label "CPS" (`name_source` member, `name_xref` `s_CPS_cyt_pa`), 16 nodes: 10 species, 6 reactions. Species: CPS, Geranylgeranyl-PP, KAO1,2, KO, KS, Phytoene, ent-Copalyl-PP, ent-Kaurene, ent-Kaurenoic acid, vsiRNA34327. The reactions of this model have no names.
- **M12** (smallest), label "SUR1" (`name_source` member, `name_xref` `s_SUR1_cyt_pa`), 4 nodes: 3 species, 1 reactions. Species: S-alkyl-thiohydroximate, SUR1, Thiohydroximate. The reactions of this model have no names.

### skm_ja

199 nodes. Degree bound of R2: 7. Candidates: 3; cut as currency species by R3: 1. The currency step is the same for the three variants.

| candidate | degree | share of its reactions stranded | result |
|---|---|---|---|
| MYC2 (`s_MYC2_nuc_pa`) | 21 | 0.67 | stays in the graph (hub) |
| JAZ (`s_JAZ_nuc_pa`) | 9 | 0.33 | currency (cut, cloned) |
| JAM (`s_JAM_cyt_pa`) | 8 | 0.50 | stays in the graph (hub) |

| | leiden | infomap | walktrap |
|---|---|---|---|
| groups L1 / L2 / L3 | 12 / 0 / 0 | 12 / 0 / 0 | 11 / 0 / 0 |
| L1 size min / median / max | 10 / 14 / 48 | 8 / 13 / 51 | 6 / 18 / 55 |
| L2 size min / median / max | - | - | - |
| L3 size min / median / max | - | - | - |
| leaf groups (no children): number | 12 | 12 | 11 |
| leaf size min / median / max | 10 / 14 / 48 | 8 / 13 / 51 | 6 / 18 / 55 |
| largest L1 group, share of all nodes | 24.1% | 25.6% | 27.6% |
| leaf groups below 5 nodes / above 60 nodes | 0 / 0 | 0 / 0 | 0 / 0 |
| unnamed groups (share) | 0 of 12 (0%) | 0 of 12 (0%) | 0 of 11 (0%) |
| name sources | member 8, notes 2, sbo 2 | member 8, notes 1, sbo 3 | member 7, notes 1, sbo 3 |
| duplicate labels | 0 | 0 | 0 |
| modules returned by the algorithm (R4, before R5-R6) | 16 | 41 | 14 |
| modularity: R4 output / L1 / leaf groups | 0.844 / 0.828 / 0.828 | 0.775 / 0.828 / 0.828 | 0.824 / 0.806 / 0.806 |
| nodes placed by the algorithm itself (rule R4 or R7) | 89% | 46% | 94% |
| nodes per rule | R4 178, R5b 10, R6 10, R8 1 | R4 91, R5a 37, R5b 10, R6 60, R8 1 | R4 188, R5b 10, R8 1 |
| currency species with clones / clones in total | 1 / 4 | 1 / 2 | 1 / 2 |
| leaf groups lying in one compartment (`loc_space`) | 8% | 17% | 27% |
| nodes in the majority compartment of their leaf group | 73% | 76% | 73% |

Agreement between the variants (adjusted Rand index / normalised mutual information over all nodes):

| pair | L1: ARI / NMI | leaf groups: ARI / NMI |
|---|---|---|
| leiden vs infomap | 0.78 / 0.84 | 0.78 / 0.84 |
| leiden vs walktrap | 0.76 / 0.83 | 0.76 / 0.83 |
| infomap vs walktrap | 0.81 / 0.86 | 0.81 / 0.86 |

Stability (seeds 1 to 10, adjusted Rand index over the 45 pairs of runs; the delivered files are seed 1):

| | L1: mean (min) | leaf groups: mean (min) | modules from R4, range | modularity of R4 output, range |
|---|---|---|---|---|
| leiden | 0.86 (0.62) | 0.86 (0.62) | 14 to 16 | 0.844 to 0.846 |
| infomap | 0.82 (0.62) | 0.82 (0.62) | 39 to 42 | 0.771 to 0.782 |
| walktrap | no random step: identical for every seed (checked, seeds 1 and 7) | same | - | - |

Walktrap instead depends on the walk length. Against the delivered 4 steps: 3 steps ARI 0.71 (L1) and 0.86 (leaf groups); 5 steps ARI 0.86 (L1) and 0.86 (leaf groups).

Three example groups of `topology_leiden` (the largest, the median and the smallest leaf group):

- **M01** (largest), label "genetic production" (`name_source` sbo, `name_xref` `SBO:0000589`), 48 nodes: 26 species, 22 reactions. Species: CLH, DELLA, DELLA|JAZ, DELLA|MYC, EIN3(like), EIN3(like)|JAZ, ICS, JAZ, JAZ|WRKY57, JR1, and others. The reactions of this model have no names.
- **M07** (median), label "active transport" (`name_source` sbo, `name_xref` `SBO:0000657`), 14 nodes: 9 species, 5 reactions. Species: CA, CTS, JASSY, OPDA, OPR, TGA. The reactions of this model have no names.
- **M12** (smallest), label "MYC2" (`name_source` member, `name_xref` `s_MYC2_cyt_pa`), 10 nodes: 6 species, 4 reactions. Species: EDR1, MPK15, MYC2, PP2AB1. The reactions of this model have no names.

### skm_sa

352 nodes. Degree bound of R2: 8. Candidates: 4; cut as currency species by R3: 2. The currency step is the same for the three variants.

| candidate | degree | share of its reactions stranded | result |
|---|---|---|---|
| SARD1 (`s_SARD1_nuc_pa`) | 31 | 0.71 | stays in the graph (hub) |
| SA (`s_SA_cyt_m`) | 16 | 0.38 | currency (cut, cloned) |
| ICS (`s_ICS_chl_pa`) | 13 | 0.15 | currency (cut, cloned) |
| NPR1 (`s_NPR1_nuc_pa`) | 10 | 0.50 | stays in the graph (hub) |

| | leiden | infomap | walktrap |
|---|---|---|---|
| groups L1 / L2 / L3 | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 8 / 2 |
| L1 size min / median / max | 19 / 28 / 48 | 20 / 28 / 48 | 5 / 20 / 90 |
| L2 size min / median / max | - | - | 6 / 11 / 64 |
| L3 size min / median / max | - | - | 22 / 32 / 42 |
| leaf groups (no children): number | 12 | 12 | 19 |
| leaf size min / median / max | 19 / 28 / 48 | 20 / 28 / 48 | 5 / 14 / 42 |
| largest L1 group, share of all nodes | 13.6% | 13.6% | 25.6% |
| leaf groups below 5 nodes / above 60 nodes | 0 / 0 | 0 / 0 | 0 / 0 |
| unnamed groups (share) | 0 of 12 (0%) | 0 of 12 (0%) | 0 of 22 (0%) |
| name sources | member 12 | member 12 | member 20, notes 1, sbo 1 |
| duplicate labels | 0 | 0 | 0 |
| modules returned by the algorithm (R4, before R5-R6) | 23 | 68 | 26 |
| modularity: R4 output / L1 / leaf groups | 0.857 / 0.840 / 0.840 | 0.769 / 0.843 / 0.843 | 0.795 / 0.780 / 0.829 |
| nodes placed by the algorithm itself (rule R4 or R7) | 76% | 23% | 56% |
| nodes per rule | R4 267, R5b 12, R6 71, R8 2 | R4 80, R5a 55, R5b 12, R6 203, R8 2 | R4 128, R5a 92, R5b 12, R6 49, R7 69, R8 2 |
| currency species with clones / clones in total | 2 / 9 | 2 / 8 | 2 / 12 |
| leaf groups lying in one compartment (`loc_space`) | 0% | 0% | 10% |
| nodes in the majority compartment of their leaf group | 64% | 66% | 72% |

Agreement between the variants (adjusted Rand index / normalised mutual information over all nodes):

| pair | L1: ARI / NMI | leaf groups: ARI / NMI |
|---|---|---|
| leiden vs infomap | 0.75 / 0.85 | 0.75 / 0.85 |
| leiden vs walktrap | 0.40 / 0.67 | 0.58 / 0.77 |
| infomap vs walktrap | 0.43 / 0.72 | 0.58 / 0.77 |

Stability (seeds 1 to 10, adjusted Rand index over the 45 pairs of runs; the delivered files are seed 1):

| | L1: mean (min) | leaf groups: mean (min) | modules from R4, range | modularity of R4 output, range |
|---|---|---|---|---|
| leiden | 0.88 (0.71) | 0.88 (0.71) | 22 to 24 | 0.856 to 0.857 |
| infomap | 0.81 (0.69) | 0.81 (0.69) | 66 to 68 | 0.768 to 0.773 |
| walktrap | no random step: identical for every seed (checked, seeds 1 and 7) | same | - | - |

Walktrap instead depends on the walk length. Against the delivered 4 steps: 3 steps ARI 0.55 (L1) and 0.88 (leaf groups); 5 steps ARI 0.48 (L1) and 0.67 (leaf groups).

Three example groups of `topology_leiden` (the largest, the median and the smallest leaf group):

- **M01** (largest), label "SARD1" (`name_source` member, `name_xref` `s_SARD1_nuc_pa`), 48 nodes: 23 species, 25 reactions. Species: ADR1, ADR1-L1, ADR1-L2, ALD1, BAK1, BAP1, BAP2, BIK1, BON1, CDPK, and others. The reactions of this model have no names.
- **M07** (median), label "SAGT" (`name_source` member, `name_xref` `s_SAGT_cyt_pa`), 27 nodes: 13 species, 14 reactions. Species: ICS, NAC019, NAC019|TCP8, NAC055, NAC072, SAG, SAGT, SARD1|TCP8, TCP8, TCP8|WRKY28, and others. The reactions of this model have no names.
- **M12** (smallest), label "EDS1" (`name_source` member, `name_xref` `s_EDS1_cyt_p`), 19 nodes: 10 species, 9 reactions. Species: CAMTA3, EDS1, EDS1-WRKY18, EIJP1, ICS, NRG1.1, NRG1.1-EDS1, TNL-effector. The reactions of this model have no names.

## 2. The three algorithms side by side

Everything except R4 is shared, so the currency species, the size bounds and the naming are identical. Modularity is measured on the graph without currency species.

| model | modules from R4: leiden / infomap / walktrap | modularity of R4 output | modularity of leaf groups | nodes placed by the algorithm itself | largest L1 share |
|---|---|---|---|---|---|
| brenda_merged | 22 / 47 / 31 | 0.795 / 0.741 / 0.756 | 0.789 / 0.761 / 0.763 | 83% / 48% / 71% | 14% / 16% / 16% |
| e_coli_core | 20 / 31 / 37 | 0.703 / 0.668 / 0.627 | 0.692 / 0.688 / 0.656 | 82% / 56% / 54% | 20% / 22% / 32% |
| R-HSA-77289 | 16 / 24 / 22 | 0.799 / 0.756 / 0.762 | 0.775 / 0.762 / 0.767 | 78% / 48% / 59% | 23% / 16% / 19% |
| skm_pss | 68 / 261 / 133 | 0.925 / 0.834 / 0.879 | 0.894 / 0.889 / 0.883 | 83% / 45% / 66% | 16% / 16% / 20% |
| skm_ja | 16 / 41 / 14 | 0.844 / 0.775 / 0.824 | 0.828 / 0.828 / 0.806 | 89% / 46% / 94% | 24% / 26% / 28% |
| skm_sa | 23 / 68 / 26 | 0.857 / 0.769 / 0.795 | 0.840 / 0.843 / 0.829 | 76% / 23% / 56% | 14% / 14% / 26% |

What the table shows:

- **Number of modules.** Leiden returns the fewest modules before repair on five models (22, 20, 16, 68, 16, 23 in panel order); Infomap returns the most on five (47, 31, 24, 261, 41, 68); Walktrap (31, 37, 22, 133, 14, 26) returns the most on e_coli_core and the fewest on skm_ja. On the directed graph most reactions are irreversible and many species are dead ends, so Infomap's flow gets trapped in small pieces.
- **Modularity.** Leiden's own output has the highest modularity on all six models. After R5 to R7 the leaf groups of Leiden still lead on five; on skm_sa Infomap is ahead by 0.003. The repair rules raise modularity above the algorithm's own output on 6 of 6 models for Infomap and 5 of 6 for Walktrap, and on 0 for Leiden: for the first two the repair, not the algorithm, did part of the work.
- **How much of the result is the algorithm's.** Share of nodes whose group was fixed by R4 or R7 and not moved by R5 or R6: Leiden 76% to 89%, Walktrap 54% to 94%, Infomap 23% to 56%. The Infomap variant is to a large part a product of the merge rules.
- **Size.** Largest L1 share: Leiden 14% to 24%, Infomap 14% to 26%, Walktrap 16% to 32%. No variant has a catch-all group. Walktrap is the least even.
- **Agreement.** The variants agree moderately (section 1): ARI of leaf groups between 0.41 and 0.81 depending on model and pair. On skm_pss the leaf groups agree well (ARI 0.69 to 0.80) while the L1 families do not (ARI 0.23 to 0.34): the small modules are robust, the way they are merged into 12 families is not.
- **Stability.** Over seeds 1 to 10 the leaf groups have a mean ARI of 0.62 to 0.88 for Leiden and 0.60 to 0.82 for Infomap. L1 on skm_pss drops to 0.42 for both. Walktrap has no random step and is identical for every seed, but it moves with its one parameter: against 4 steps, 3 or 5 steps give a leaf-group ARI of 0.48 to 0.93.

## 3. Recommendation

**`topology_leiden`.** It needs the least repair, so the delivered groups are closest to what the algorithm found; it has the highest modularity before repair on all six models and after repair on five; its largest top-level group is the smallest of the three, or tied, on five of six models (R-HSA-77289 is the exception); and it returns a number of modules closest to what a page holds. Its weakness is the seed: two seeds give visibly different groupings on the metabolic models (leaf-group ARI near 0.63 on brenda_merged and R-HSA-77289). The delivered result is one of several near-equal optima, fixed by seed 1.

`topology_walktrap` is the choice when a grouping must not depend on a seed at all; it pays with less even groups and a lower modularity. `topology_infomap` is not recommended on these models: its modules are too small to stand, and the result owes more to R5 and R6 than to flow.

## 4. Known failure cases

1. **Most names are hub names, not functions.** Over the six models, 149 of 167 Leiden groups are named after a member (Infomap 156 of 176, Walktrap 165 of 186). An annotation term qualified under R9a for only 18 Leiden groups (ec 1, go 1, notes 14, sbo 2). A module found by connectivity rarely coincides with one EC class, GO term, SBO term or MapMan bin. A hub name says what the group is centred on, not what it does. No group is unnamed and no label occurs twice.
2. **GO and EC labels are bare codes.** The feature tables carry no GO or EC names, so those groups are labelled `GO:0006635` or `EC 6.2.1`. The identifier resolves, but a reader cannot read it.
3. **The top level of the large model mixes pathways.** skm_pss has 68 Leiden modules; R6 merges them into 8 families of 179 to 295 nodes, each named after one hub: VPg (295), CAX1 (258), ORA59 (233), GA (221), CBP60G (209), NPR1 (203), IAA (191), cZ (179). The label describes one part of such a family. The L2 groups inside are the readable level.
4. **Small separate components take top-level places.** In skm_pss 4 of the 12 L1 groups are separate components of 12, 7, 5 and 4 nodes next to the families of 179 to 295. In R-HSA-77289 two components of 4 nodes are L1 groups. They cannot be merged on connectivity (R5c), and they stay below S_MIN.
5. **Compartments are ignored.** Share of Leiden leaf groups that lie in one compartment: brenda_merged 100%, e_coli_core 18%, R-HSA-77289 100%, skm_pss 36%, skm_ja 8%, skm_sa 0%. In the models with several compartments most groups straddle a compartment border, and a drawing that nests compartments has to split them. (100% for the single-compartment brenda_merged and the nearly single-compartment Reactome pathway says nothing.)
6. **Residual glue in E. coli core.** The degree bound is 11. NADP+, NADPH, CO2 and coenzyme A have degree 9 and stay in the graph. One Leiden group is centred on NADP+ and holds glutamate metabolism together with other NADP-dependent steps.
7. **Main metabolites cut as currency.** By brenda_merged's own reaction note `main_metabolites`, 12 of the 15 species cut there are side metabolites in most of their reactions; L-glutamate, AMP and 2-oxoglutarate are main metabolites and are cut all the same. In E. coli core pyruvate is cut. These species get one home and clones, which is tolerable for junction metabolites but is not what 'currency' means. This check exists for one model only.
8. **The glue test depends on order and on context.** In skm_pss SARD1 (31 reactions) is cut and CBP60G (30) stays, although the two regulate the same targets: once SARD1 is gone, the targets hang on CBP60G alone. MYC2 stays as a hub in skm_ja (0.67 of its reactions stranded) and is cut in skm_pss (0.24), where its targets have other connections. In the signalling models the species that are cut are regulators and messengers (Ca2+, SA, ROS), not currency in the metabolic sense.
9. **A star cannot be split.** The largest leaf group of skm_pss has 133 nodes (bound: 60): the viral protein VPg with its 53 reactions. Its family already carries the label VPg, so the star itself got the next hub name, "6K2|VPg", which is misleading. One Walktrap leaf group of brenda_merged has 61 nodes because currency species are added after the size check.
10. **The page bound overrides the algorithm.** In all 18 runs the algorithm returned more than 12 modules. In 15 of them parts of at least S_MIN nodes had to be merged by R6 to reach 12; in the other 3 the merging of small parts (R5) was enough. Either way the top level is not the algorithm's own partition.
11. **Seed dependence** of Leiden and Infomap, and step-length dependence of Walktrap, as in section 2. No consensus over seeds was built.
12. **All links weigh 1.** Stoichiometry, the difference between a substrate and a modifier, and the sign of a regulation play no part.

## 5. Limits of these notes

- The rules were written while looking at these six models. No rule has a per-model branch or number, but the choice of the Poisson bound over Tukey's fences and the definition of 'stranded' were made after seeing all six. They are untested on other models.
- Walktrap's nested levels come from running Walktrap again inside a large group (R7), the same way as for Leiden and Infomap, not from cutting its global dendrogram at several heights.
- Modularity before and after the currency cut is not comparable: removing nodes raises it by itself. In an exploratory Leiden run on the uncut graph (seed 1, not delivered) modularity was 0.568 (brenda_merged), 0.520 (e_coli_core), 0.581 (R-HSA-77289), 0.889 (skm_pss), 0.825 (skm_ja), 0.820 (skm_sa).
- No curated labels were used or seen. Whether the groups match pathways is not measured here.
