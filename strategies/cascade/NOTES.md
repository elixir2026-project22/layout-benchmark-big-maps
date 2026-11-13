# cascade - notes

Results of `strategy.py` on the six panel models, produced from the feature tables of the bundle only.
No reference labels were read; nothing below is a score. All counts were read back from the files in `out/`.

## Checks

`python common/contract.py validate` and two runs of `strategy.py` (PYTHONHASHSEED 0 and 98765), compared by
SHA-256 of `groups.tsv` and `groups_meta.tsv`:

| model | validate | run 1 vs run 2 | groups.tsv sha256 (first 16) | groups_meta.tsv sha256 (first 16) |
|---|---|---|---|---|
| brenda_merged | valid | identical | `8d12d91bc5bdd007` | `8b18e9e3f3b63e91` |
| e_coli_core | valid | identical | `b6a6657f9b22c944` | `1ceb9e20d703222f` |
| R-HSA-77289 | valid | identical | `57792b5254ef625d` | `63a7fd64c11657d1` |
| skm_pss | valid | identical | `6c8d2ba4ff14912a` | `0841cb02f6d6a34d` |
| skm_ja | valid | identical | `798bd550602c6ef9` | `643fc6ae272d5664` |
| skm_sa | valid | identical | `ea1c69de066dd186` | `fe4a49e4168a7374` |

`out/<model>/<model>_grouped.xml` was written with `common/contract.py write-sbml --strategy cascade` for all
six models. libSBML reads R-HSA-77289_grouped.xml with one error (6020305, a `<layout>` without `layout:id`);
the same error is in the input file `panel/models/R-HSA-77289.xml` and does not come from the groups.

Run time: under 0.5 s per model.

## Per model

Levels: level 1 = anchor family, route, currency module, or unconnected steps; level 2 = stage or route
element (or the localisation level, where a group has no stage or element); level 3 = localisation.
Group ids: `C..` anchor family, `.s1/.s2/.s3` its stages (production, molecular interaction, genetic
production), `P..` route, `.c` cycle, `.b1` linear branch, `K01` currency module, `U_<space>` unconnected
steps, `<group>.<compartment id>` localisation level.

## brenda_merged

396 nodes (178 species, 218 reactions). CAP = 15. Anchor families: 0. Currency species: 19.

| level | groups | nodes covered | size min / median / max | groups under 5 nodes |
|---|---|---|---|---|
| 1 | 22 | 396 (100%) | 4 / 17 / 40 | 2 |
| 2 | 13 | 163 (41%) | 5 / 10 / 26 | 0 |
| 3 | 0 | 0 | - | - |

- Top-level sizes: 40, 38, 37, 29, 26, 24, 24, 22, 22, 19, 18, 16, 13, 13, 12, 9, 8, 7, 6, 5, 4, 4. Largest top-level group: 10.1% of the nodes.
- Unnamed groups: 0 of 35 (0%). Name sources: compartment 1, member 34. Bases: flow 25, localisation 1, topology 9.
- Nodes per rule: R10 174, R11 39, R13 3, R14 2, R17b 19, R17c 127, R17d 32.
- Species with clones: 51.
- Top-level groups (id, label, nodes): P01 GDP / succinyl-CoA (40); P02 alpha-D-ribose 1-phosphate / adenosine (38); P03 5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxylate / 5-aminoimidazole-4-carboxamide ribonucleotide (37); P04 D-fructose 6-phosphate / D-fructose 1,6-bisphosphate (29); P05 6-phosphogluconate / D-glucose (26); P06 dGTP / 2-deoxy-D-ribose 1-phosphate (24); P07 (S)-malate / isocitrate (24); K01 H2O / ATP / ADP (22); P08 D-glycerate / glutaconyl-CoA (22); P09 5-amino-4-imidazole carboxylate / L-serine (19); P10 Carbamoyl phosphate / L-citrulline (18); P11 3-phosphoglycerate / 2-phospho-D-glycerate (16); P13 (R)-4'-phosphopantothenate / 4'-phosphopantetheine (13); P12 succinate semialdehyde / 4-aminobutanoate (13); P14 6,8-dihydroxypurine / purine (12); P15 7-methylxanthine / Caffeine (9); P16 (3S)-citryl-CoA / citrate (8); P17 (+)-citramalate / mesaconate (7); P18 D-ribulose 1,5-bisphosphate / alpha-D-ribose 1,5-bisphosphate (6); P19 acetyl phosphate / glycine (5); U_c default (4); P20 dATP (4).

Three example groups:

- `P08` **D-glycerate / glutaconyl-CoA** (level 1, basis flow, named from member `M_mt_996ef3e8d7`): 22 nodes = 12 species + 10 reactions. Species: glutaconyl-CoA, Butyrate, acetyl-CoA, crotonyl-CoA, (R)-2-hydroxyglutarate, D-glycerate, oxidized acceptor, (R)-2-hydroxyglutaryl-1-CoA, reduced acceptor, butyryl-CoA, acceptor, acetate. Reactions: BR50094, BS370908, BS358809, NODE_10018, NODE_11298, NODE_2957, NODE_8287, BR51640 ….
- `P10` **Carbamoyl phosphate / L-citrulline** (level 1, basis flow, named from member `M_mt_73617a0182`): 18 nodes = 8 species + 10 reactions. Species: allophanate, Carbamoyl phosphate, HCO3-, L-citrulline, argininosuccinate, L-ornithine, Urea, L-arginine. Reactions: BR142, BR47825, BR19238, BR3431, BR19260, BR43950, BR143, BR22698 ….
- `P12` **succinate semialdehyde / 4-aminobutanoate** (level 1, basis flow, named from member `M_mt_c20619d98d`): 13 nodes = 5 species + 8 reactions. Species: L-alanine, 4-hydroxybutyrate, 4-aminobutanoate, succinate semialdehyde, 2-Oxoglutaramate. Reactions: BR51069, BR29388, BR47437, BR47673, BR275, BR50321, BR29388-a, BR50304.

## e_coli_core

167 nodes (72 species, 95 reactions). CAP = 10. Anchor families: 0. Currency species: 20.

| level | groups | nodes covered | size min / median / max | groups under 5 nodes |
|---|---|---|---|---|
| 1 | 12 | 167 (100%) | 4 / 7 / 35 | 1 |
| 2 | 6 | 99 (59%) | 8 / 16 / 25 | 0 |
| 3 | 0 | 0 | - | - |

- Top-level sizes: 35, 32, 28, 22, 10, 7, 7, 6, 6, 5, 5, 4. Largest top-level group: 21.0% of the nodes.
- Unnamed groups: 0 of 18 (0%). Name sources: compartment 4, member 14. Bases: flow 12, localisation 4, topology 2.
- Nodes per rule: R10 53, R11 27, R13 15, R17b 20, R17c 36, R17d 16.
- Species with clones: 27.
- Top-level groups (id, label, nodes): K01 H+ / H2O H2O / ATP C10H12N5O13P3 (35); P01 CO2 CO2 / L-Malate (32); P02 3-Phospho-D-glycerate / D-Erythrose 4-phosphate (28); P03 Ammonium / L-Glutamate (22); P04 Acetaldehyde / Ethanol (10); P05 Acetate / Acetyl phosphate (7); P06 Fumarate (7); P08 D-Fructose 1,6-bisphosphate / Dihydroxyacetone phosphate (6); P07 Formate (6); P10 6-Phospho-D-gluconate / 6-phospho-D-glucono-1,5-lactone (5); P09 D-Lactate (5); P11 Succinate (4).

Three example groups:

- `P01.c` **CO2 CO2 / L-Malate / cycle through Isocitrate** (level 2, basis topology, named from member `M_icit_c`): 19 nodes = 8 species + 11 reactions. Species: Cis-Aconitate, 2-Oxoglutarate, Citrate, CO2 CO2, Glyoxylate, Isocitrate, Oxaloacetate, Succinate. Reactions: Aconitase (half-reaction A, C…, Aconitase (half-reaction B, I…, Citrate synthase, Isocitrate dehydrogenase (NAD…, Isocitrate lyase, Malate synthase, Malate dehydrogenase, Malic enzyme (NAD) ….
- `P03` **Ammonium / L-Glutamate** (level 1, basis flow, named from member `M_nh4_c`): 22 nodes = 8 species + 14 reactions. Species: 2-Oxoglutarate, L-Glutamine, L-Glutamine, L-Glutamate, L-Glutamate, Ammonium, Ammonium, Succinyl-CoA. Reactions: 2-Oxogluterate dehydrogenase, 2 oxoglutarate reversible tra…, 2-Oxoglutarate exchange, L-Glutamine exchange, L-Glutamate exchange, Ammonia exchange, Glutamine synthetase, L-glutamine transport via ABC… ….
- `K01` **H+ / H2O H2O / ATP C10H12N5O13P3** (level 1, basis flow, named from member `M_h_c`): 35 nodes = 20 species + 15 reactions. Species: Acetyl-CoA, ADP C10H12N5O10P2, ATP C10H12N5O13P3, Coenzyme A, H2O H2O, H2O H2O, H+, H+, Nicotinamide adenine dinucleotide, Nicotinamide adenine dinucleotide - r…, Nicotinamide adenine dinucleotide pho…, Nicotinamide adenine dinucleotide pho… …. Reactions: ATP maintenance requirement, ATP synthase (four protons fo…, Cytochrome oxidase bd (ubiqui…, H2O exchange, H+ exchange, O2 exchange, Phosphate exchange, Pyruvate exchange ….

## R-HSA-77289

148 nodes (101 species, 47 reactions). CAP = 7. Anchor families: 0. Currency species: 8.

| level | groups | nodes covered | size min / median / max | groups under 5 nodes |
|---|---|---|---|---|
| 1 | 11 | 148 (100%) | 4 / 15 / 24 | 1 |
| 2 | 0 | 0 | - | - |
| 3 | 0 | 0 | - | - |

- Top-level sizes: 24, 21, 17, 17, 15, 15, 12, 8, 8, 7, 4. Largest top-level group: 16.2% of the nodes.
- Unnamed groups: 0 of 11 (0%). Name sources: compartment 3, member 8. Bases: flow 8, localisation 3.
- Nodes per rule: R10 37, R11 2, R14 8, R17b 8, R17c 51, R17d 42.
- Species with clones: 18.
- Top-level groups (id, label, nodes): P01 ATP [mitochondrial matrix] / AMP [mitochondrial matrix] (24); U_compartment_5460 mitochondrial matrix (21); P02 4-cis-decenoyl-CoA [mitochondrial matrix] / 2-trans-4-cis-decadienoyl-CoA [mitochondrial matrix] (17); P03 DEC-CoA [mitochondrial matrix] / (S)-3-Hydroxydodecanoyl-CoA [mitochondrial matrix] (17); P04 MYS-CoA [mitochondrial matrix] / trans-Tetradec-2-enoyl-CoA [mitochondrial matrix] (15); P05 3-Oxohexanoyl-CoA [mitochondrial matrix] / (S)-Hydroxyhexanoyl-CoA [mitochondrial matrix] (15); U_compartment_70101 cytosol (12); K01 H+ [mitochondrial matrix] / CoA-SH [mitochondrial matrix] / FAD [mitochondrial matrix] (8); P06 Octanoyl-CoA [mitochondrial matrix] / 3-Oxodecanoyl-CoA [mitochondrial matrix] (8); P07 3HB-CoA [mitochondrial matrix] / Crotonoyl-CoA [mitochondrial matrix] (7); U_compartment_17957 endoplasmic reticulum lumen (4).

Three example groups:

- `P04` **MYS-CoA [mitochondrial matrix] / trans-Tetradec-2-enoyl-CoA [mitochondrial matrix]** (level 1, basis flow, named from member `species_77272`): 15 nodes = 9 species + 6 reactions. Species: Trifunctional Protein [mitochondrial …, MYS-CoA [mitochondrial matrix], trans-Tetradec-2-enoyl-CoA [mitochond…, (S)-3-Hydroxytetradecanoyl-CoA [mitoc…, ACADVL dimer [mitochondrial matrix], PalmCoA [mitochondrial matrix], trans-Hexadec-2-enoyl-CoA [mitochondr…, (S)-3-Hydroxyhexadecanoyl-CoA [mitoch…, 3-Oxopalmitoyl-CoA [mitochondrial mat…. Reactions: myristoyl-CoA+FAD => trans-te…, trans-Tetradec-2-enoyl-CoA+H2…, palmitoyl-CoA+FAD => trans-He…, trans-Hexadec-2-enoyl-CoA+H2O…, (S)-3-Hydroxyhexadecanoyl-CoA…, 3-Oxopalmitoyl-CoA+CoA-SH => ….
- `P02` **4-cis-decenoyl-CoA [mitochondrial matrix] / 2-trans-4-cis-decadienoyl-CoA [mitochondrial matrix]** (level 1, basis flow, named from member `species_109986`): 17 nodes = 11 species + 6 reactions. Species: 4-cis-decenoyl-CoA [mitochondrial mat…, 2-trans-4-cis-decadienoyl-CoA [mitoch…, 3-trans-decenoyl-CoA [mitochondrial m…, DCI dimer [mitochondrial matrix], DECR1 tetramer [mitochondrial matrix], ACADL tetramer [mitochondrial matrix], ACADM tetramer [mitochondrial matrix], tdec2-CoA [mitochondrial matrix], cis,cis-3,6-Dodecadienoyl-CoA [mitoch…, LIN-CoA [mitochondrial matrix], trans,cis-Lauro-2,6-dienoyl-CoA  [mit…. Reactions: Isomerization of cis,cis-3,6-…, Linoleoyl-CoA is converted to…, dehydrogenation of 4-cis-dece…, Trans,cis-lauro-2,6-dienoyl-C…, Reduction of 2-trans-4-cis-de…, Isomerization of 3-trans-dece….
- `U_compartment_5460` **mitochondrial matrix** (level 1, basis localisation, named from compartment `compartment_5460`): 21 nodes = 16 species + 5 reactions. Species: acyl-CoA [mitochondrial matrix], NDUFAB1 [mitochondrial inner membrane], ACOT2,9, THEM4,5 dimer [mitochondrial…, MCFA, LCFA [mitochondrial matrix], MCFA-CoA, LCFA-CoA [mitochondrial mat…, S-2MPDA-CoA [mitochondrial matrix], S-2MDPDA-CoA [mitochondrial matrix], ACAD10 [mitochondrial matrix], BH-CoA [mitochondrial matrix], DHB-CoA [mitochondrial matrix], ACAD11 [mitochondrial inner membrane], ACAA2 tetramer [mitochondrial matrix] …. Reactions: ACOT2,9,THEM4,5 hydrolyse MCF…, ACAD10 dehydrogenates S-2MPDA…, ACAD11 dehydrogenates BH-CoA, ACAA2 tetramer transfers acyl…, MCAT transfers Mal from Mal-C….

## skm_pss

1817 nodes (961 species, 856 reactions). CAP = 30. Anchor families: 21. Currency species: 0.

| level | groups | nodes covered | size min / median / max | groups under 5 nodes |
|---|---|---|---|---|
| 1 | 24 | 1817 (100%) | 7 / 52 / 281 | 0 |
| 2 | 41 | 1711 (94%) | 6 / 30 / 229 | 0 |
| 3 | 54 | 1349 (74%) | 5 / 13 / 220 | 0 |

- Top-level sizes: 281, 240, 218, 178, 123, 96, 86, 72, 67, 59, 57, 53, 52, 51, 41, 38, 28, 16, 14, 12, 10, 10, 8, 7. Largest top-level group: 15.5% of the nodes.
- Unnamed groups: 0 of 119 (0%). Name sources: compartment 61, member 23, sbo 35. Bases: flow 58, localisation 61.
- Nodes per rule: R5 156, R6 77, R7 559, R8 26, R10 17, R12 16, R14 5, R17a 26, R17c 618, R17d 317.
- Species with clones: 132.
- Top-level groups (id, label, nodes): C01 Ca2+ (281); C04 potyvirus (240); C13 ABA (218); C02 SA (178); C16 JA-Ile (123); C05 DZ / cZ / iP / +1 (96); C06 ROS (86); C09 Cu2+ / ET (72); C15 IAA (67); C08 red light (59); C03 Heat (57); C17 OG (53); C10 GA (52); C14 Brassinolide (51); C11 MeCLA / SL (41); C20 p-Coumaric acid (38); P01 CA-CoA / BA-CoA (28); U_c_cyt cytoplasm (16); C07 Waterlogging (14); P02 InsP5 / InsP6 (12); C19 eATP (10); C21 trichous-bacteria (10); C12 cholesterol (8); C18 UDP-Glc (7).

Three example groups:

- `C16.s1.c_per` **JA-Ile / production / peroxisome** (level 3, basis localisation, named from compartment `c_per`): 26 nodes = 16 species + 10 reactions. Species: ACH, ACX, ACX, CTS, JA-CoA, JAT2, JA, KAT, MFP, OPC4-CoA, OPC6-CoA, OPC8-CoA …. Reactions: rx00031, rx00033, rx00034, rx00035, rx00036, rx00037, rx00038, rx00039 ….
- `C05.s2` **DZ / cZ / iP / +1 / molecular interaction** (level 2, basis flow, named from sbo `SBO:0000344`): 26 nodes = 12 species + 14 reactions. Species: AHK2,3,4, AHK2,3,4, AHP, AHP, AHP, ARR-A, ARR-A, ARR-B, ARR-B, MYB2, MYB2, PR1. Reactions: rx00384, rx00385, rx00386, rx00387, rx00388, rx00389, rx00390, rx00391 ….
- `C09` **Cu2+ / ET** (level 1, basis flow, named from member `s_Cu2_er_m`): 72 nodes = 38 species + 34 reactions. Species: ABI4, ABI4, ACC, ACO, ACO, ACS, ACS, ACS, CTR|ETR, CTR|ETR, CTR, Cu2+ …. Reactions: activation_SAHH_pa_er, activation_SAMS_pa_er, rx00001, rx00002, rx00003, rx00004, rx00005, rx00006 ….

## skm_ja

199 nodes (117 species, 82 reactions). CAP = 10. Anchor families: 1. Currency species: 0.

| level | groups | nodes covered | size min / median / max | groups under 5 nodes |
|---|---|---|---|---|
| 1 | 4 | 199 (100%) | 10 / 54 / 81 | 0 |
| 2 | 9 | 196 (98%) | 5 / 26 / 38 | 0 |
| 3 | 0 | 0 | - | - |

- Top-level sizes: 81, 67, 41, 10. Largest top-level group: 40.7% of the nodes.
- Unnamed groups: 0 of 13 (0%). Name sources: compartment 9, member 1, sbo 3. Bases: flow 4, localisation 9.
- Nodes per rule: R5 19, R6 5, R7 34, R8 12, R10 3, R12 9, R17a 1, R17c 73, R17d 43.
- Species with clones: 7.
- Top-level groups (id, label, nodes): C01.s1 JA-Ile: production (81); C01.s2 JA-Ile: molecular interaction (67); C01.s3 JA-Ile: genetic production (41); P01 SAHH / L-homo-cys (10).

Three example groups:

- `C01.s1.c_chl` **JA-Ile: production / chloroplast** (level 2, basis localisation, named from compartment `c_chl`): 18 nodes = 11 species + 7 reactions. Species: 12,13-EOT, 13-HPOT, ALA, AOC, AOS, DAD1, LOX, MDA, OPDA, ROS, phospholipids. Reactions: rx00014, rx00026, rx00029, rx00030, rx00396, rx00693, transport_LOX_pa_cyt_to_chl.
- `C01.s2.c_nuc` **JA-Ile: molecular interaction / nucleus** (level 2, basis localisation, named from compartment `c_nuc`): 29 nodes = 12 species + 17 reactions. Species: COI1, DELLA|JAZ, DELLA|MYC, EIN3(like)|JAZ, EIN3(like), JAM, MYC2, MYC2, MYC|PYL, WD/bHLH/MYB, WD/bHLH/MYB, WRKY57. Reactions: activation_WDbHLHMYB_ca_nuc, rx00018, rx00043, rx00044, rx00047, rx00062, rx00128, rx00159 ….
- `C01.s3` **JA-Ile: genetic production** (level 1, basis flow, named from sbo `SBO:0000589`): 41 nodes = 22 species + 19 reactions. Species: CLH, CML|Ca2+, COI1, DELLA, ERF1, ICS, JAZ, JR1, LOX2|VPg, MKK4,5, NAC019, NAC055 …. Reactions: rx00046, rx00049, rx00050, rx00051, rx00052, rx00163, rx00196, rx00206 ….

## skm_sa

352 nodes (193 species, 159 reactions). CAP = 13. Anchor families: 4. Currency species: 0.

| level | groups | nodes covered | size min / median / max | groups under 5 nodes |
|---|---|---|---|---|
| 1 | 7 | 352 (100%) | 6 / 38 / 99 | 0 |
| 2 | 11 | 320 (91%) | 7 / 22 / 78 | 0 |
| 3 | 0 | 0 | - | - |

- Top-level sizes: 99, 93, 66, 38, 37, 13, 6. Largest top-level group: 28.1% of the nodes.
- Unnamed groups: 0 of 18 (0%). Name sources: compartment 11, member 4, sbo 3. Bases: flow 7, localisation 11.
- Nodes per rule: R5 21, R6 20, R7 53, R8 17, R10 13, R12 35, R17a 4, R17c 120, R17d 69.
- Species with clones: 26.
- Top-level groups (id, label, nodes): C01.s2 SA: molecular interaction (99); C01.s3 SA: genetic production (93); C01.s1 SA: production (66); C04 p-Coumaric acid (38); P01 CA-CoA / BA-CoA (37); C02 ROS (13); C03 cholesterol (6).

Three example groups:

- `C01.s1.c_chl` **SA: production / chloroplast** (level 2, basis localisation, named from compartment `c_chl`): 19 nodes = 4 species + 15 reactions. Species: Chor, EDS5, ICS, IsoChor. Reactions: rx00069, rx00187, rx00188, rx00189, rx00190, rx00240, rx00254, rx00255 ….
- `C04` **p-Coumaric acid** (level 1, basis flow, named from member `s_pCoumaricacid_cyt_m`): 38 nodes = 24 species + 14 reactions. Species: 4CL, 4CL, ADT, CA, CA, CM, CYP73A5, DOPAL, ELI5, L-Dopa, L-arogenate, L-arogenate …. Reactions: rx00073, rx00076, rx00460, rx00540, rx00541, rx00542, rx00543, rx00695 ….
- `P01` **CA-CoA / BA-CoA** (level 1, basis flow, named from member `s_CACoA_per_m`): 37 nodes = 22 species + 15 reactions. Species: 3H3PP-CoA, 3O3PP-CoA, AAO, ABA, ADH, Abscisic aldehyde, BA-CoA, BA-CoA, BA, BA, BD, BZO1 …. Reactions: rx00036, rx00530, rx00531, rx00532, rx00533, rx00534, rx00535, rx00536 ….

## Known failure cases

Stated as they are in the output files.

1. **Route names are member species, not pathway names.** A route is labelled with its two most
   characteristic species ("Carbamoyl phosphate / L-citrulline" for the urea-cycle route of brenda_merged).
   The label says where the box is, not what the pathway is called. In e_coli_core the route holding the
   citrate/glyoxylate cycle is labelled "CO2 CO2 / L-Malate": CO2 takes part in 9 reactions, CAP is 10, and it has
   no carrier partner, so R2 does not find it and it is the busiest species of that route.
2. **Currency detection is approximate.** Against the `main_metabolites` notes of brenda_merged (not used by
   the strategy, compared afterwards): of the reaction participants that the notes do not list as main,
   347 are currency by R2 and 107 are not (recall 0.76); of the participants R2 calls currency, 32
   are listed as main (precision 0.92). Currency set of brenda_merged: ADP, ATP, CO2, H+, H2O, L-glutamate, L-glutamine, NAD+, NADH, NADP+, NADPH, NH3, O2, formaldehyde, oxidized ferredoxin, phosphate, reduced ferredoxin, thioredoxin, thioredoxin disulfide.
   L-glutamate/L-glutamine and formaldehyde are in it although they are main metabolites of some reactions.
   In e_coli_core pyruvate takes part in 11 reactions, one more than CAP, and is therefore currency; its uptake
   and exchange (PYRt2, EX_pyr_e) sit in the currency module. O2 becomes currency by R2c (in CYTBD, a
   transport-class step, it is the only non-currency participant), so CYTBD, O2t and EX_o2_e sit there too.
   Currency module of e_coli_core: ATPM, ATPS4r, CYTBD, EX_h2o_e, EX_h_e, EX_o2_e, EX_pi_e, EX_pyr_e, H2Ot, NADH16, NADTRHD, O2t, PIt2r, PYRt2, THD2.
3. **Hubs cut textbook pathways.** Routes end at branch points that join many reactions. In e_coli_core
   glycolysis is in three groups: PGI and GAPD in "3-Phospho-D-glycerate / D-Erythrose 4-phosphate", PFK and FBA in "D-Fructose 1,6-bisphosphate / Dihydroxyacetone phosphate",
   PYK in "CO2 CO2 / L-Malate". Seven of its twelve top-level groups hold seven nodes or fewer:
   Acetate / Acetyl phosphate (4 reactions), Fumarate (5 reactions), Formate (4 reactions), D-Fructose 1,6-bisphosphate / Dihydroxyacetone phosphate (4 reactions), D-Lactate (3 reactions), 6-Phospho-D-gluconate / 6-phospho-D-glucono-1,5-lactone (3 reactions), Succinate (3 reactions).
4. **Cycles are only as good as the reversibility flags.** R16 lets a reaction with `reversible` = 1 be passed
   either way. In brenda_merged 213 of 218 reactions are reversible, so any closed loop of three species
   counts: 9 of 20 routes get a cycle element and it holds 77 of their 125 reactions. Two
   parallel routes between the same metabolites are reported as a cycle. Route elements cover 41 % of the
   nodes of brenda_merged; the other nodes sit directly in their route.
5. **Unconnected steps get a location, not a meaning.** Reactions tied to the rest of the model through
   currency species only, or through nothing, are grouped by `loc_space` (R14): 8 reactions in
   R-HSA-77289 (endoplasmic reticulum lumen (4 nodes), mitochondrial matrix (21 nodes), cytosol (12 nodes)), 2 in brenda_merged, 5 in skm_pss, none in the other three. Such a group is named
   after its compartment and stands beside flow groups that lie in the same compartment.
6. **Reactions without SBO terms.** All 47 reactions of R-HSA-77289 lack an SBO term; R1 classes them from
   their participants (5 with a protein or complex among reactants or products are class state, the others
   metabolic). No anchor is found there, so a Reactome signalling pathway in this encoding would be grouped
   by the route fallback of R10, not by cascades. R-HSA-77289 has no level 2: no route meets the conditions of R16
   (a cycle and a rest of three reactions each) or of R20 (two spaces with five nodes each).
7. **One home per reaction under crosstalk.** A step reachable from several anchors goes to the single family
   with the largest flow share (R7); the other families show nothing for it, only shared species are cloned
   (R18; 132 species in skm_pss). The share is computed with an even split at every branching, which is a
   convention, not a measured flux.
8. **Stage labels are coarse.** "production" (SBO:0000393) also holds turnover and export of the anchor and
   the supply of its enzymes, including their genetic production. "genetic production" holds every step after
   the first genetic reaction on the best path: in skm_pss the potyvirus group `C04.s3` has 229 nodes and
   66 of its 87 reactions are binding reactions, not genetic productions.
9. **Top level of skm_pss has 24 groups**, more than a page shows comfortably. The largest is Ca2+
   (281 nodes, 15.5%); the smallest are UDP-Glc (7), cholesterol (8), eATP (10), trichous-bacteria (10), InsP5 / InsP6 (12), Waterlogging (14), cytoplasm (16).
   Families of three to six reactions are kept as anchors by R9.
10. **R19 is a threshold.** skm_ja and skm_sa each have one family above half of the nodes, so its stages become
    the top level ("JA-Ile: production", ...). A model whose largest family holds 49 % would keep it as one box
    with the stages at level 2.
11. **Small top-level groups remain.** Groups under 5 nodes at level 1: brenda_merged 2, e_coli_core 1,
    R-HSA-77289 1, the SKM models 0. A route of three reactions can keep as few as four nodes when its species are
    made elsewhere.
12. **Production cores are not bounded.** R5 follows mass flow into an anchor as far as it goes. In the SKM
    models this is the biosynthesis chain; in a model that combined whole metabolism with signalling it would
    pull every precursor route into the anchor's family.
13. **`--sbml` is not opened.** Everything the rules read is in the feature tables; the argument exists for
    the common CLI.
14. **Not checked against curated labels.** By design of the task no reference grouping was available to
    this track; agreement with curated pathways is unknown.
