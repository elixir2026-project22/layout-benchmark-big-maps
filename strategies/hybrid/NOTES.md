# NOTES - strategy `hybrid`

All numbers below were computed from the files in `out/<model>/` (`groups.tsv`, `groups_meta.tsv`) after the final run; sizes count species and reactions, clones are not counted. "Leaf" = the deepest group of a node. Environment of the run: Python 3.11.16, python-igraph 1.0.0, leidenalg 0.12.0, Windows.

## Overview: do the bounds hold?

| model | nodes | L1 groups (aim) | largest L1 share (cap 1/3) | largest leaf (cap 40) | groups L1 / L2 / L3 | leaves | leaves < 3 nodes | singleton leaves | unnamed groups | validate | two runs identical |
|---|---|---|---|---|---|---|---|---|---|---|---|
| brenda_merged | 396 | 10 (10) | 0.139 | 40 | 10 / 6 / 0 | 13 | 0 | 0 | 0 of 16 | valid | yes |
| e_coli_core | 167 | 5 (5) | 0.269 | 35 | 5 / 14 / 2 | 15 | 0 | 0 | 0 of 21 | valid | yes |
| R-HSA-77289 | 148 | 4 (4) | 0.331 | 39 | 4 / 5 / 2 | 8 | 0 | 0 | 0 of 11 | valid | yes |
| skm_pss | 1817 | 12 (12) | 0.201 | 40 | 12 / 59 / 115 | 138 | 21 | 12 | 0 of 186 | valid | yes |
| skm_ja | 199 | 5 (5) | 0.327 | 35 | 5 / 17 / 0 | 17 | 3 | 2 | 0 of 22 | valid | yes |
| skm_sa | 352 | 9 (9) | 0.295 | 36 | 9 / 25 / 19 | 38 | 7 | 4 | 0 of 53 | valid | yes |

Bounds stated in RULES.md: 4 to 12 level-1 groups, no level-1 group above one third of the model, no leaf group above 40 nodes. All three hold on all six models. The lower size bound (3 nodes) does not hold for localisation groups, which no rule may merge (see "Known failure cases", item 5).

## How each level was used

| model | L1 function (notes term) | L1 connectivity (member name) | function groups below L1 | reaction-type groups | localisation groups | of which same-compartment parts (R13) | connectivity parts below L1 | nodes placed by R3 / R4 / R5 / R6 / R7 | hub species cloned (copies) |
|---|---|---|---|---|---|---|---|---|---|
| brenda_merged | 0 | 10 | 0 | 0 | 0 | 0 | 6 | 0 / 0 / 278 / 0 / 118 | 19 (124) |
| e_coli_core | 0 | 5 | 0 | 3 | 13 | 0 | 0 | 0 / 0 / 75 / 0 / 92 | 4 (16) |
| R-HSA-77289 | 0 | 4 | 0 | 0 | 5 | 0 | 2 | 0 / 0 / 56 / 0 / 92 | 13 (40) |
| skm_pss | 8 | 4 | 23 | 9 | 120 | 9 | 22 | 316 / 1085 / 63 / 0 / 353 | 32 (142) |
| skm_ja | 5 | 0 | 0 | 0 | 17 | 0 | 0 | 33 / 71 / 22 / 0 / 73 | 4 (5) |
| skm_sa | 6 | 3 | 4 | 3 | 36 | 0 | 1 | 59 / 156 / 62 / 0 / 75 | 13 (34) |

## brenda_merged

- **Nodes:** 396 (178 species, 218 reactions). **Groups:** 16 (level 1: 10, level 2: 6, level 3: 0). Nodes whose deepest group is at level 1 / 2 / 3: 242 / 154 / 0.
- **Groups by level and basis:** L1 topology: 10, L2 topology: 6.
- **Size distribution, level 1:** 55, 51, 48, 40, 39, 39, 35, 33, 32, 24 nodes; largest share 0.139.
- **Size distribution, leaf groups (13):** min 21, lower quartile 26, median 28, upper quartile 35, max 40; 0 below 3 nodes, of which 0 singletons. All groups: median 33, max 55; most children under one group: 2.
- **Unnamed groups:** 0 of 16 (0 %). Name sources: member 16.
- **Level 1 construction:** 0 nodes carry function terms, 0 of them confirmed as seeds (R3); 23 connectivity communities (R2); 23 groups before merging, 13 merges (R7) down to 10.
- **Groups above 40 nodes without finer function terms (R10/R11 order):** M01 `oxaloacetate`: neither a type nor a localisation split was admissible; connectivity parts (R12); M02 `alpha-D-ribose 1-phosphate`: neither a type nor a localisation split was admissible; connectivity parts (R12); M03 `D-glucose 6-phosphate`: neither a type nor a localisation split was admissible; connectivity parts (R12).

| level-1 group | label | named from | nodes | share | label fit: F1 / precision / recall (carriers in group) |
|---|---|---|---|---|---|
| M01 | oxaloacetate | member: `M_mt_e7c775ae46` | 55 | 0.139 | n/a (member name) |
| M02 | alpha-D-ribose 1-phosphate | member: `M_mt_8edfd22915` | 51 | 0.129 | n/a (member name) |
| M03 | D-glucose 6-phosphate | member: `M_mt_c40006a428` | 48 | 0.121 | n/a (member name) |
| M04 | 3-phosphoglycerate | member: `M_mt_1137269b50` | 40 | 0.101 | n/a (member name) |
| M05 | thioredoxin | member: `M_mt_d1fa415e8a` | 39 | 0.098 | n/a (member name) |
| M06 | Carbamoyl phosphate | member: `M_mt_73617a0182` | 39 | 0.098 | n/a (member name) |
| M07 | 5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxylate | member: `M_mt_1a547d42f7` | 35 | 0.088 | n/a (member name) |
| M08 | D-fructose 6-phosphate | member: `M_mt_dfef3a54ea` | 33 | 0.083 | n/a (member name) |
| M09 | succinate semialdehyde | member: `M_mt_c20619d98d` | 32 | 0.081 | n/a (member name) |
| M10 | glutaconyl-CoA | member: `M_mt_12580d2e8b` | 24 | 0.061 | n/a (member name) |

**Three example groups**

- *largest level-1 group* `M01` - "oxaloacetate" (level 1, basis topology, 55 nodes, 2 child groups): 24 species: `succinate`, `oxaloacetate`, `acetyl-CoA`, `fumarate`, `CoA`, `(S)-malate`, `succinyl-CoA`, `citrate`, `GDP`, `GTP` (+14 more); 31 reactions: `BR47671`, `BR53166`, `BR3376`, `BR50544` (+27 more).
- *median-size leaf group* `M01.p1` - "oxaloacetate / succinyl-CoA" (level 2, basis topology, 28 nodes): 13 species: `succinyl-CoA`, `succinate`, `GDP`, `fumarate`, `GTP`, `L-aspartate`, `cGMP`, `adenylosuccinate`, `Maleate`, `FADH2` (+3 more); 15 reactions: `BR47671`, `BR53166`, `BR50544`, `BR88` (+11 more).
- *smallest leaf group* `M03.p2` - "D-glucose 6-phosphate / 5-amino-4-imidazole carboxylate" (level 2, basis topology, 21 nodes): 11 species: `5-amino-4-imidazole carboxylate`, `L-serine`, `5,10-methenyltetrahydrofolate`, `5-formimidoyltetrahydrofolate`, `tetrahydrofolate`, `5,10-methylenetetrahydrofolate`, `glycine`, `Formiminoglycine`, `imidazolinone`, `5-ureido-4-imidazole carboxylate` (+1 more); 10 reactions: `BR40573`, `BR36902`, `BR43943`, `BR3426` (+6 more).

## e_coli_core

- **Nodes:** 167 (72 species, 95 reactions). **Groups:** 21 (level 1: 5, level 2: 14, level 3: 2). Nodes whose deepest group is at level 1 / 2 / 3: 0 / 153 / 14.
- **Groups by level and basis:** L1 topology: 5, L2 localisation: 11, L2 reaction_type: 3, L3 localisation: 2.
- **Size distribution, level 1:** 45, 39, 37, 24, 22 nodes; largest share 0.269.
- **Size distribution, leaf groups (15):** min 3, lower quartile 6, median 8, upper quartile 17, max 35; 0 below 3 nodes, of which 0 singletons. All groups: median 12, max 45; most children under one group: 3.
- **Unnamed groups:** 0 of 21 (0 %). Name sources: compartment 13, member 5, sbo 3.
- **Level 1 construction:** 0 nodes carry function terms, 0 of them confirmed as seeds (R3); 20 connectivity communities (R2); 20 groups before merging, 15 merges (R7) down to 5.
- **Groups above 40 nodes without finer function terms (R10/R11 order):** M01 `Ammonium`: type first (leaf groups below 3 nodes: 0 with type first, 0 with localisation first); children by reaction_type.

| level-1 group | label | named from | nodes | share | label fit: F1 / precision / recall (carriers in group) |
|---|---|---|---|---|---|
| M01 | Ammonium | member: `M_nh4_c` | 45 | 0.269 | n/a (member name) |
| M02 | D-Fructose 6-phosphate | member: `M_f6p_c` | 39 | 0.234 | n/a (member name) |
| M03 | Acetaldehyde | member: `M_acald_c` | 37 | 0.222 | n/a (member name) |
| M04 | Formate | member: `M_for_c` | 24 | 0.144 | n/a (member name) |
| M05 | Fumarate | member: `M_fum_c` | 22 | 0.132 | n/a (member name) |

**Three example groups**

- *largest level-1 group* `M01` - "Ammonium" (level 1, basis topology, 45 nodes, 3 child groups): 20 species: `H+`, `ATP C10H12N5O13P3`, `H2O H2O`, `ADP C10H12N5O10P2`, `Phosphate`, `L-Glutamate`, `L-Glutamine`, `Nicotinamide adenine dinucleotide phosphate`, `Nicotinamide adenine dinucleotide phosphate - reduced`, `Ammonium` (+10 more); 25 reactions: `Adenylate kinase`, `2 oxoglutarate reversible transport via symport`, `ATP maintenance requirement`, `ATP synthase (four protons for one ATP)` (+21 more).
- *median-size leaf group* `M01.sbo0000185.c_b` - "Ammonium / translocation reaction / cytosol (boundary)" (level 3, basis localisation, 8 nodes): no species; 8 reactions: `2 oxoglutarate reversible transport via symport`, `ATP synthase (four protons for one ATP)`, `L-glutamine transport via ABC system`, `L glutamate transport via proton symport  reversible` (+4 more).
- *smallest leaf group* `M04.c_b` - "Formate / cytosol (boundary)" (level 2, basis localisation, 3 nodes): no species; 3 reactions: `Acetate reversible transport via proton symport`, `Formate transport via diffusion`, `Formate transport in via proton symport`.

## R-HSA-77289

- **Nodes:** 148 (101 species, 47 reactions). **Groups:** 11 (level 1: 4, level 2: 5, level 3: 2). Nodes whose deepest group is at level 1 / 2 / 3: 67 / 35 / 46.
- **Groups by level and basis:** L1 topology: 4, L2 localisation: 5, L3 topology: 2.
- **Size distribution, level 1:** 49, 39, 32, 28 nodes; largest share 0.331.
- **Size distribution, leaf groups (8):** min 3, lower quartile 12, median 17, upper quartile 28, max 39; 0 below 3 nodes, of which 0 singletons. All groups: median 28, max 49; most children under one group: 3.
- **Unnamed groups:** 0 of 11 (0 %). Name sources: compartment 5, member 6.
- **Level 1 construction:** 0 nodes carry function terms, 0 of them confirmed as seeds (R3); 17 connectivity communities (R2); 17 groups before merging, 13 merges (R7) down to 4.
- **Groups above 40 nodes without finer function terms (R10/R11 order):** M01 `DEC-CoA [mitochondrial matrix]`: type first (leaf groups below 3 nodes: 0 with type first, 0 with localisation first); children by localisation.

| level-1 group | label | named from | nodes | share | label fit: F1 / precision / recall (carriers in group) |
|---|---|---|---|---|---|
| M01 | DEC-CoA [mitochondrial matrix] | member: `species_77308` | 49 | 0.331 | n/a (member name) |
| M02 | ATP [mitochondrial matrix] | member: `species_113593` | 39 | 0.264 | n/a (member name) |
| M03 | MYS-CoA [mitochondrial matrix] | member: `species_77272` | 32 | 0.216 | n/a (member name) |
| M04 | 4-cis-decenoyl-CoA [mitochondrial matrix] | member: `species_109986` | 28 | 0.189 | n/a (member name) |

**Three example groups**

- *largest level-1 group* `M01` - "DEC-CoA [mitochondrial matrix]" (level 1, basis topology, 49 nodes, 2 child groups): 32 species: `CoA-SH [mitochondrial matrix]`, `Ac-CoA [mitochondrial matrix]`, `H2O [mitochondrial matrix]`, `ECHS1 hexamer [mitochondrial matrix]`, `HADH dimer [mitochondrial matrix]`, `DEC-CoA [mitochondrial matrix]`, `tdec2-CoA [mitochondrial matrix]`, `(S)-3-Hydroxydodecanoyl-CoA [mitochondrial matrix]`, `3-Oxododecanoyl-CoA [mitochondrial matrix]`, `2-trans-Dodecenoyl-CoA [mitochondrial matrix]` (+22 more); 17 reactions: `ACOT2,9,THEM4,5 hydrolyse MCFA-CoA, LCFA-CoA`, `(S)-3-Hydroxydodecanoyl-CoA+NAD => 3-Oxododecanoyl-CoA+NADH+H`, `2-trans-Dodecenoyl-CoA+H2O =>( S)-3-Hydroxydodecanoyl-CoA`, `lauroyl-CoA+FAD => 2-trans-Dodecenoyl-CoA+FADH2` (+13 more).
- *median-size leaf group* `M01.compartment_5460.p2` - "DEC-CoA [mitochondrial matrix] / mitochondrial matrix / 3-Oxooctanoyl-CoA [mitochondrial matrix]" (level 3, basis topology, 17 nodes): 13 species: `CoA-SH [mitochondrial matrix]`, `Ac-CoA [mitochondrial matrix]`, `3-Oxooctanoyl-CoA [mitochondrial matrix]`, `acyl-CoA [mitochondrial matrix]`, `ACOT2,9, THEM4,5 dimer [mitochondrial matrix]`, `MCFA, LCFA [mitochondrial matrix]`, `MCFA-CoA, LCFA-CoA [mitochondrial matrix]`, `HX-CoA [mitochondrial matrix]`, `(S)-Hydroxyoctanoyl-CoA [mitochondrial matrix]`, `ACAA2 tetramer [mitochondrial matrix]` (+3 more); 4 reactions: `ACOT2,9,THEM4,5 hydrolyse MCFA-CoA, LCFA-CoA`, `3-Oxooctanoyl-CoA+CoA-SH => Hexanoyl-CoA + acetyl-CoA`, `(S)-Hydroxyoctanoyl-CoA+NAD => 3-Oxooctanoyl-CoA+NADH+H`, `ACAA2 tetramer transfers acyl group from Ac-CoA to acyl-CoA forming 3OA-CoA and CoA-SH`.
- *smallest leaf group* `M01.compartment_5460_b` - "DEC-CoA [mitochondrial matrix] / mitochondrial matrix (boundary)" (level 2, basis localisation, 3 nodes): 2 species: `NDUFAB1 [mitochondrial inner membrane]`, `Mal-PPANT-S112-NDUFAB1 [mitochondrial inner membrane]`; 1 reactions: `MCAT transfers Mal from Mal-CoA to NDUFAB1`.

## skm_pss

- **Nodes:** 1817 (961 species, 856 reactions). **Groups:** 186 (level 1: 12, level 2: 59, level 3: 115). Nodes whose deepest group is at level 1 / 2 / 3: 0 / 344 / 1473.
- **Groups by level and basis:** L1 function: 8, L1 topology: 4, L2 function: 23, L2 localisation: 17, L2 reaction_type: 4, L2 topology: 15, L3 localisation: 103, L3 reaction_type: 5, L3 topology: 7.
- **Size distribution, level 1:** 365, 350, 187, 151, 149, 115, 103, 94, 93, 86, 83, 41 nodes; largest share 0.201.
- **Size distribution, leaf groups (138):** min 1, lower quartile 3, median 9, upper quartile 22, max 40; 21 below 3 nodes, of which 12 singletons (all of basis localisation). All groups: median 17, max 365; most children under one group: 10.
- **Unnamed groups:** 0 of 186 (0 %). Name sources: compartment 120, member 26, notes 31, sbo 9.
- **Level 1 construction:** 528 nodes carry function terms, 316 of them confirmed as seeds (R3); 110 connectivity communities (R2); 75 groups before merging, 63 merges (R7) down to 12.
- **Groups above 40 nodes without finer function terms (R10/R11 order):** M04 `VPg`: type first (leaf groups below 3 nodes: 1 with type first, 1 with localisation first); children by reaction_type; M07 `HC-Pro`: localisation first (leaf groups below 3 nodes: 5 with type first, 0 with localisation first); children by localisation; M08 `ARF`: localisation first (leaf groups below 3 nodes: 7 with type first, 0 with localisation first); children by localisation; M11 `ROS`: localisation first (leaf groups below 3 nodes: 9 with type first, 0 with localisation first); children by localisation; M12 `Phytohormone action: gibberellin: biosynthesis`: type first (leaf groups below 3 nodes: 0 with type first, 0 with localisation first); children by reaction_type.

| level-1 group | label | named from | nodes | share | label fit: F1 / precision / recall (carriers in group) |
|---|---|---|---|---|---|
| M01 | Phytohormone action | notes: `mapman:10` | 365 | 0.201 | 0.69 / 0.78 / 0.63 (79) |
| M02 | RNA biosynthesis | notes: `mapman:14` | 350 | 0.193 | 0.69 / 0.74 / 0.64 (86) |
| M03 | External stimuli response | notes: `mapman:26` | 187 | 0.103 | 0.71 / 0.82 / 0.64 (54) |
| M04 | VPg | member: `s_VPg_cyt_p` | 151 | 0.083 | n/a (member name) |
| M05 | Multi-process regulation | notes: `mapman:27` | 149 | 0.082 | 0.64 / 0.80 / 0.54 (36) |
| M06 | Solute transport | notes: `mapman:24` | 115 | 0.063 | 0.67 / 0.69 / 0.64 (18) |
| M07 | HC-Pro | member: `s_HCPro_cyt_pa` | 103 | 0.057 | n/a (member name) |
| M08 | ARF | member: `s_ARF_nuc_pa` | 94 | 0.052 | n/a (member name) |
| M09 | Protein modification | notes: `mapman:18` | 93 | 0.051 | 0.68 / 0.90 / 0.55 (28) |
| M10 | Enzyme classification | notes: `mapman:50` | 86 | 0.047 | 0.22 / 0.85 / 0.13 (22) |
| M11 | ROS | member: `s_ROS_cyt_m` | 83 | 0.046 | n/a (member name) |
| M12 | Phytohormone action: gibberellin: biosynthesis | notes: `mapman:10.6.1` | 41 | 0.023 | 0.80 / 0.67 / 1.00 (6) |

**Three example groups**

- *largest level-1 group* `M01` - "Phytohormone action" (level 1, basis function, 365 nodes, 10 child groups): 190 species: `MYC2`, `IAA`, `SA`, `JAM`, `cZ`, `iP`, `tZ`, `AHK2,3,4`, `DZ`, `EIN3(like)` (+180 more); 175 reactions: `activation_JAM_pa_cyt`, `rx00002`, `rx00003`, `rx00004` (+171 more).
- *median-size leaf group* `M08.c_nuc` - "ARF / nucleus" (level 2, basis localisation, 9 nodes): 4 species: `ARF`, `ARF`, `ABI4`, `GH3.7`; 5 reactions: `activation_ARF_pa_nuc`, `rx00112`, `rx00465`, `rx00518` (+1 more).
- *smallest leaf group* `M01.f1.c_ga` - "Phytohormone action / ethylene / golgi apparatus" (level 3, basis localisation, 1 nodes): 1 species: `RTE1`; 0 reactions.

## skm_ja

- **Nodes:** 199 (117 species, 82 reactions). **Groups:** 22 (level 1: 5, level 2: 17, level 3: 0). Nodes whose deepest group is at level 1 / 2 / 3: 0 / 199 / 0.
- **Groups by level and basis:** L1 function: 5, L2 localisation: 17.
- **Size distribution, level 1:** 65, 40, 37, 33, 24 nodes; largest share 0.327.
- **Size distribution, leaf groups (17):** min 1, lower quartile 4, median 10, upper quartile 17, max 35; 3 below 3 nodes, of which 2 singletons (all of basis localisation). All groups: median 12, max 65; most children under one group: 4.
- **Unnamed groups:** 0 of 22 (0 %). Name sources: compartment 17, notes 5.
- **Level 1 construction:** 65 nodes carry function terms, 33 of them confirmed as seeds (R3); 34 connectivity communities (R2); 21 groups before merging, 16 merges (R7) down to 5.
- **Groups above 40 nodes without finer function terms (R10/R11 order):** M01 `RNA biosynthesis`: localisation first (leaf groups below 3 nodes: 4 with type first, 0 with localisation first); children by localisation.

| level-1 group | label | named from | nodes | share | label fit: F1 / precision / recall (carriers in group) |
|---|---|---|---|---|---|
| M01 | RNA biosynthesis | notes: `mapman:14` | 65 | 0.327 | 0.72 / 0.65 / 0.81 (17) |
| M02 | Enzyme classification | notes: `mapman:50` | 40 | 0.201 | 0.45 / 0.70 / 0.33 (7) |
| M03 | Phytohormone action | notes: `mapman:10` | 37 | 0.186 | 0.42 / 0.83 / 0.28 (10) |
| M04 | Phytohormone action: jasmonic acid: biosynthesis | notes: `mapman:10.7.1` | 33 | 0.166 | 0.70 / 0.78 / 0.64 (7) |
| M05 | Enzyme classification: EC_1 oxidoreductases | notes: `mapman:50.1` | 24 | 0.121 | 0.62 / 0.62 / 0.62 (5) |

**Three example groups**

- *largest level-1 group* `M01` - "RNA biosynthesis" (level 1, basis function, 65 nodes, 3 child groups): 34 species: `MYC2`, `JAM`, `JAZ`, `JAZ`, `MYC2`, `DELLA`, `JAZ`, `JAM`, `JAM`, `MYC2` (+24 more); 31 reactions: `activation_JAM_pa_cyt`, `activation_WDbHLHMYB_ca_nuc`, `rx00018`, `rx00046` (+27 more).
- *median-size leaf group* `M02.c_per` - "Enzyme classification / peroxisome" (level 2, basis localisation, 10 nodes): 6 species: `JA-CoA`, `JA`, `KAT`, `3O3PP-CoA`, `ACH`, `JAT2`; 4 reactions: `rx00037`, `rx00038`, `rx00039`, `rx00535`.
- *smallest leaf group* `M04.c_nuc` - "Phytohormone action: jasmonic acid: biosynthesis / nucleus" (level 2, basis localisation, 1 nodes): 1 species: `TGA`; 0 reactions.

## skm_sa

- **Nodes:** 352 (193 species, 159 reactions). **Groups:** 53 (level 1: 9, level 2: 25, level 3: 19). Nodes whose deepest group is at level 1 / 2 / 3: 0 / 149 / 203.
- **Groups by level and basis:** L1 function: 6, L1 topology: 3, L2 function: 4, L2 localisation: 20, L2 topology: 1, L3 localisation: 16, L3 reaction_type: 3.
- **Size distribution, level 1:** 104, 70, 58, 28, 24, 22, 16, 16, 14 nodes; largest share 0.295.
- **Size distribution, leaf groups (38):** min 1, lower quartile 3, median 8, upper quartile 14, max 36; 7 below 3 nodes, of which 4 singletons (all of basis localisation). All groups: median 12, max 104; most children under one group: 5.
- **Unnamed groups:** 0 of 53 (0 %). Name sources: compartment 36, member 4, notes 10, sbo 3.
- **Level 1 construction:** 114 nodes carry function terms, 59 of them confirmed as seeds (R3); 47 connectivity communities (R2); 26 groups before merging, 17 merges (R7) down to 9.
- **Groups above 40 nodes without finer function terms (R10/R11 order):** M02 `Phytohormone action`: localisation first (leaf groups below 3 nodes: 4 with type first, 2 with localisation first); children by localisation.

| level-1 group | label | named from | nodes | share | label fit: F1 / precision / recall (carriers in group) |
|---|---|---|---|---|---|
| M01 | External stimuli response | notes: `mapman:26` | 104 | 0.295 | 0.66 / 0.55 / 0.83 (24) |
| M02 | Phytohormone action | notes: `mapman:10` | 70 | 0.199 | 0.63 / 0.64 / 0.61 (16) |
| M03 | RNA biosynthesis | notes: `mapman:14` | 58 | 0.165 | 0.58 / 0.76 / 0.47 (16) |
| M04 | Enzyme classification | notes: `mapman:50` | 28 | 0.080 | 0.22 / 0.67 / 0.13 (4) |
| M05 | EDS5 | member: `s_EDS5_cyt_pa` | 24 | 0.068 | n/a (member name) |
| M06 | Secondary metabolism: phenolics biosynthesis: p-coumaroyl-CoA biosynthesis | notes: `mapman:8.2.1` | 22 | 0.062 | 0.80 / 0.67 / 1.00 (4) |
| M07 | 3H3PP-CoA | member: `s_3H3PPCoA_per_m` | 16 | 0.045 | n/a (member name) |
| M08 | Amino acid metabolism | notes: `mapman:4` | 16 | 0.045 | 0.67 / 1.00 / 0.50 (3) |
| M09 | AAO | member: `s_AAO_cyt_pa` | 14 | 0.040 | n/a (member name) |

**Three example groups**

- *largest level-1 group* `M01` - "External stimuli response" (level 1, basis function, 104 nodes, 3 child groups): 52 species: `SARD1`, `CBP60G`, `EDS1`, `NPR1`, `EDS1`, `EDS1`, `PAD4`, `PAD4`, `SARD1`, `CAT` (+42 more); 52 reactions: `activation_CAT_pa_per`, `activation_EDS1_pa_cyt`, `activation_GST_pa_cyt`, `rx00087` (+48 more).
- *median-size leaf group* `M06.c_per` - "Secondary metabolism: phenolics biosynthesis: p-coumaroyl-CoA biosynthesis / peroxisome" (level 2, basis localisation, 8 nodes): 4 species: `CA`, `4CL`, `CA-CoA`, `CTS`; 4 reactions: `rx00531`, `rx00532`, `rx00533`, `transport_4CL_pa_cyt_to_per`.
- *smallest leaf group* `M01.f2.c_chl` - "External stimuli response / pathogen: defense mechanisms: systemic acquired resistance (SAR) / chloroplast" (level 3, basis localisation, 1 nodes): no species; 1 reactions: `rx00190`.

## Known failure cases

Stated plainly; the numbers are from the final run.

1. **The reaction-type level is used rarely.** The brief puts reaction type at level 2 everywhere. Here it is applied only to a group above 40 nodes that has no finer function terms, and only where it does not shatter the localisation level. Result: reaction-type groups per model: brenda_merged 0, e_coli_core 3, R-HSA-77289 0, skm_pss 9, skm_ja 0, skm_sa 3. `brenda_merged` has one SBO term for all reactions and `R-HSA-77289` has no SBO terms on reactions, so the level cannot exist there at all; in `skm_ja` every module was either within the leaf bound or produced more tiny groups with type first. A reader looking for "activation / binding / genetic production" boxes will mostly not find them.
2. **In three of six models level 1 is pure connectivity, named after one member.** `brenda_merged`, `e_coli_core` and `R-HSA-77289` carry no function terms that the SBML also names (EC numbers are reaction types; the GO ids on Reactome reactions come without names and are mostly molecular functions), so every level-1 group is a merged Leiden community labelled with one metabolite (R14). Such a label says where the module is, not what it does: `ATP [mitochondrial matrix]` (R-HSA-77289 M02) holds propionyl-CoA catabolism, the C4 and C6 rounds of beta-oxidation and two CoA ligations, `Ammonium` (e_coli_core M01) mixes glutamate metabolism with ATP maintenance and PEP carboxylation. Level-1 groups named after a member: brenda_merged 10 of 10, e_coli_core 5 of 5, R-HSA-77289 4 of 4, skm_pss 4 of 12, skm_ja 0 of 5, skm_sa 3 of 9.
3. **Generic MapMan bins name modules.** A term names a group when a majority of the annotated members carry it. The bin `50 Enzyme classification` therefore names 4 level-1 groups (skm_pss M10: F1 0.22, recall 0.13; skm_ja M02: F1 0.45, recall 0.33; skm_ja M05: F1 0.62, recall 0.62; skm_sa M04: F1 0.22, recall 0.13): true for the members, but it describes an enzyme class, not a pathway, and most carriers of the term sit elsewhere. `14 RNA biosynthesis` groups are transcription-factor families collected from all over the network, not one pathway. In a single-pathway model the pathway's own term is the background that R3 subtracts, so the core of `skm_ja` (MYC2, JAZ, JAM and their targets, 65 nodes) is labelled `RNA biosynthesis`, not jasmonic acid signalling.
4. **Function modules are not connected pieces of the network.** Seeds keep their term wherever they are (R3), so a level-1 function group is a union of islands. Level-1 groups that are not internally connected: brenda_merged 0 of 10, e_coli_core 0 of 5, R-HSA-77289 1 of 4, skm_pss 10 of 12, skm_ja 3 of 5, skm_sa 3 of 9. Worst case `skm_pss` M06 `Solute transport`: 6 components, the largest holds 35 of 115 nodes. In `R-HSA-77289` the model itself has 4 disconnected components and the three small ones were merged into M03 without any connecting edge (R7, third case).
5. **Tiny localisation groups remain.** Localisation is given and a localisation group cannot be merged with another compartment, so a module with one or two nodes in an organelle gets a group of one or two nodes. Leaf groups below 3 nodes: brenda_merged 0 of 13 (0 singletons), e_coli_core 0 of 15 (0 singletons), R-HSA-77289 0 of 8 (0 singletons), skm_pss 21 of 138 (12 singletons), skm_ja 3 of 17 (2 singletons), skm_sa 7 of 38 (4 singletons). All of them are localisation groups. This is far below the current pipeline (167 singletons of 346 groups on the SKM model) but not zero.
6. **Reaction-only and species-only boxes.** A boundary group may hold only boundary nodes, which are transport reactions; their substrates sit in the neighbouring space groups. Leaf groups holding only reactions / only species: brenda_merged 0 / 0, e_coli_core 5 / 1, R-HSA-77289 0 / 0, skm_pss 25 / 6, skm_ja 3 / 2, skm_sa 6 / 2.
7. **The same compartment several times under one parent (R13).** Where function, then localisation, still leaves more than 40 nodes and no level is left, the compartment is cut into connectivity parts labelled `cytoplasm [member]`. This happens for 9 groups, all in `skm_pss`. The parts are valid localisation groups but a drawing shows two or more "cytoplasm" boxes side by side inside one module.
8. **Merging puts unrelated material under one label.** R7 folds small groups into larger ones until the level-1 aim is met. Share of nodes whose level-1 group was decided by a merge: brenda_merged 30%, e_coli_core 55%, R-HSA-77289 62%, skm_pss 19%, skm_ja 37%, skm_sa 21%. The label is re-evaluated after merging (R8) and only says what the majority of annotated members are; the label-fit columns above show how well.
9. **Hormone pathways are level 2, not level 1, in `skm_pss`, and salicylic acid has no group of its own.** With at most 12 top-level groups the MapMan top bins win; ethylene, jasmonic acid, auxin, strigolactone, gibberellin and cytokinin are children of M01 `Phytohormone action`. Three of its ten children are connectivity groups named after a member (`tZ-ribotide`, `BB`, `BZR2`). The 10 species carrying `10.8 salicylic acid` end up in 6 different level-2 groups under 4 level-1 modules (M01, M02, M03, M06), so no group is named after salicylic acid.
10. **Only one annotation source seeds modules.** The notes key `mapman` is the only hierarchical, named function vocabulary in the panel. KEGG pathway cross-references on SKM metabolites, GO ids and EC numbers are not used for level 1. A model whose function annotation sits only in such identifiers is grouped as if it had none.
11. **Currency metabolites are down-weighted, not removed.** Edge weight 1/degree lets ATP, NAD, H2O still pull a little; `e_coli_core` M01 is the visible result. Hubs are cloned (R15), which hides the cross-links in a drawing but not in the grouping: share of species-reaction edges crossing level-1 borders: brenda_merged 36% (without hub species 10%), e_coli_core 31% (without hub species 23%), R-HSA-77289 30% (without hub species 3%), skm_pss 11% (without hub species 6%), skm_ja 8% (without hub species 2%), skm_sa 12% (without hub species 2%).
12. **Bounds that are not guaranteed.** The lower bound of 4 level-1 groups is an aim: a model with fewer than four connectivity communities (for instance one of 40 nodes or fewer) gets fewer groups; no such model is in the panel. The 3-node minimum has the exception of item 5. The number of children under one group is not bounded by a rule; the observed maximum is 10 (`skm_pss`, children of M01).
13. **`rule` in groups.tsv records level 1 only.** It holds the rule that decided the node's level-1 module (R3, R4, R5, R6 or R7). The deeper placement follows from R9-R13 and can be read from the `basis` of the node's groups. R6 was never needed at level 1 on this panel (no group exceeded a third of a model before merging); it fired inside modules at level 2 (term descent, e.g. `14.5`, `26.11`, `27.8`).
14. **Determinism rests on Leiden with a fixed seed.** Two runs, also with different `PYTHONHASHSEED`, give byte-identical files on this machine. A different leidenalg or igraph version may give different communities and therefore different groups; this was not tested.
15. **Not compared with any reference.** This track had no curated labels. Whether the modules agree with curated pathways is unknown; every statement above about "sensible" modules rests on reading member names.
