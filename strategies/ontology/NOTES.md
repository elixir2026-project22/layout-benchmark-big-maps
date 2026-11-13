# Strategy `ontology` - notes

All numbers below were computed from the `out/<model>/groups.tsv` and `groups_meta.tsv` in this folder (parameters `M_MIN` = 3, `PAGE` = 12, `S_MAX` = 0.5, `CLONE_MIN_BOXES` = 3). Sizes are in nodes (species + reactions). No curated reference label was read at any point; nothing outside the bundle was read.

## Overview

| model | species / reactions | groups L1 / L2 / L3 | largest L1 share | unnamed groups | groups < 3 nodes | reactions placed by R8 or R10 | cloned species |
|---|---|---|---|---|---|---|---|
| brenda_merged | 178 / 218 | 7 / 24 / 0 | 0.27 | 0 of 31 | 0 | 0 (0%) | 44 (25%) |
| e_coli_core | 72 / 95 | 8 / 13 / 0 | 0.23 | 0 of 21 | 0 | 1 (1%) | 25 (35%) |
| R-HSA-77289 | 101 / 47 | 6 / 6 / 0 | 0.32 | 0 of 12 | 0 | 2 (4%) | 6 (6%) |
| skm_pss | 961 / 856 | 12 / 34 / 52 | 0.26 | 0 of 98 | 0 | 159 (19%) | 29 (3%) |
| skm_ja | 117 / 82 | 9 / 11 / 5 | 0.33 | 0 of 25 | 0 | 11 (13%) | 2 (2%) |
| skm_sa | 193 / 159 | 8 / 8 / 14 | 0.32 | 0 of 30 | 0 | 25 (16%) | 6 (3%) |

## brenda_merged

- Nodes: 396 (178 species, 218 reactions). Groups: 31 (L1 7, L2 24, L3 0).
- Size distribution (nodes per group): L1: n=7, min 11, median 54, max 108, groups under 5 nodes 0; L2: n=24, min 3, median 13, max 39, groups under 5 nodes 1.
- Unnamed groups: 0 of 31 (0%). Name sources: ec 30, sbo 1. Labels that are a bare identifier (GO id, pathway id, EC number without a name): 24.
- Largest level-1 box: 0.27 of all nodes.
- Reactions by rule: R1 214, R9 4. Species by rule: R12 178.
- Localisation: 0% of nodes sit in a localisation box (R14); 100% of nodes sit in an innermost box whose members all share one `loc_space`.
- Boxes dropped as too small (R5) at level 1: 0.
- Level-1 boxes (nodes, label, identifier): 108 `EC 1 (oxidoreductases)` [ec-code:1]; 107 `EC 2 (transferases)` [ec-code:2]; 55 `EC 3 (hydrolases)` [ec-code:3]; 54 `EC 4 (lyases)` [ec-code:4]; 31 `EC 6 (ligases)` [ec-code:6]; 30 `EC 5 (isomerases)` [ec-code:5]; 11 `biochemical reaction` [SBO:0000176].

Example groups:

1. **EC 1 (oxidoreductases)** - level 1, basis function, identifier `ec-code:1`, 108 nodes (50 species, 58 reactions). Species: (S)-malate, 2-oxoglutarate, 4-hydroxybutyrate, 5,10-methylenetetrahydrofolate, 6,8-dihydroxypurine, 6-phospho-D-glucono-1,5-lactone, 6-phosphogluconate, 7-methylxanthine .... Reactions: BR50094; BS287796; BR53166; BR50577-a ....
2. **EC 4 (lyases) / EC 4.2.1** - level 2, basis function, identifier `ec-code:4.2.1`, 13 nodes (5 species, 8 reactions). Species: (R)-2-hydroxyglutaryl-1-CoA, cis-aconitate, fumarate, mesaconate, phosphoenolpyruvate. Reactions: BR3385; BR419; BR414; BR418 ....
3. **biochemical reaction** - level 1, basis reaction_type, identifier `SBO:0000176`, 11 nodes (7 species, 4 reactions). Species: 4-aminoimidazole, 5-amino-4-imidazole carboxylate, 5-ureido-4-imidazole carboxylate, 8-hydroxypurine, oxidized electron donor, purine, reduced electron acceptor. Reactions: NODE_4393; NODE_4439; NODE_4449; NODE_4455.

## e_coli_core

- Nodes: 167 (72 species, 95 reactions). Groups: 21 (L1 8, L2 13, L3 0).
- Size distribution (nodes per group): L1: n=8, min 8, median 25, max 39, groups under 5 nodes 0; L2: n=13, min 3, median 7, max 21, groups under 5 nodes 3.
- Unnamed groups: 0 of 21 (0%). Name sources: compartment 2, ec 16, sbo 3. Labels that are a bare identifier (GO id, pathway id, EC number without a name): 11.
- Largest level-1 box: 0.23 of all nodes.
- Reactions by rule: R1 49, R10 1, R9 45. Species by rule: R12 72.
- Localisation: 17% of nodes sit in a localisation box (R14); 100% of nodes sit in an innermost box whose members all share one `loc_space`.
- Boxes dropped as too small (R5) at level 1: 5.
- Level-1 boxes (nodes, label, identifier): 39 `exchange reaction` [SBO:0000627]; 30 `translocation reaction` [SBO:0000185]; 29 `EC 1 (oxidoreductases)` [ec-code:1]; 25 `EC 2 (transferases)` [ec-code:2]; 16 `EC 4 (lyases)` [ec-code:4]; 11 `EC 5 (isomerases)` [ec-code:5]; 9 `biochemical reaction` [SBO:0000176]; 8 `EC 3 (hydrolases)` [ec-code:3].

Example groups:

1. **exchange reaction** - level 1, basis reaction_type, identifier `SBO:0000627`, 39 nodes (19 species, 20 reactions). Species: 2-Oxoglutarate, Acetaldehyde, Acetate, Ammonium, CO2 CO2, D-Fructose, D-Glucose, D-Lactate .... Reactions: Acetate exchange; Acetaldehyde exchange; 2-Oxoglutarate exchange; CO2 exchange ....
2. **EC 4 (lyases) / EC 4.2.1** - level 2, basis function, identifier `ec-code:4.2.1`, 7 nodes (3 species, 4 reactions). Species: Cis-Aconitate, Citrate, H2O H2O. Reactions: Aconitase (half-reaction A, Citrate hydro-lyase); Aconitase (half-reaction B, Isocitrate hydro-lyase); Enolase; Fumarase.
3. **translocation reaction / cytosol boundary** - level 2, basis localisation, identifier `c`, 21 nodes (0 species, 21 reactions). No species has its home here. Reactions: Acetaldehyde reversible transport; Acetate reversible transport via proton symport; 2 oxoglutarate reversible transport via symport; CO2 transporter via diffusion ....

## R-HSA-77289

- Nodes: 148 (101 species, 47 reactions). Groups: 12 (L1 6, L2 6, L3 0).
- Size distribution (nodes per group): L1: n=6, min 13, median 22, max 48, groups under 5 nodes 0; L2: n=6, min 4, median 12, max 35, groups under 5 nodes 2.
- Unnamed groups: 0 of 12 (0%). Name sources: compartment 3, ec 8, go 1. Labels that are a bare identifier (GO id, pathway id, EC number without a name): 4.
- Largest level-1 box: 0.32 of all nodes.
- Reactions by rule: R1 45, R10 2. Species by rule: R12 101.
- Localisation: 14% of nodes sit in a localisation box (R14); 100% of nodes sit in an innermost box whose members all share one `loc_space`.
- Single-child chains followed (R5): `ec-code:2`; `ec-code:2.3`; `ec-code:4`; `ec-code:4.2`; `ec-code:4.2.1`.
- Boxes dropped as too small (R5) at level 1: 2.
- Level-1 boxes (nodes, label, identifier): 48 `EC 1 (oxidoreductases)` [ec-code:1]; 29 `EC 2.3.1 (transferases)` [ec-code:2.3.1]; 22 `EC 4.2.1.17 (lyases)` [ec-code:4.2.1.17]; 20 `GO:0006637` [GO:0006637]; 16 `EC 6 (ligases)` [ec-code:6]; 13 `EC 5 (isomerases)` [ec-code:5].

Example groups:

1. **EC 1 (oxidoreductases)** - level 1, basis function, identifier `ec-code:1`, 48 nodes (27 species, 21 reactions). Species: 2-trans-4-cis-decadienoyl-CoA [mitochondrial matrix], 4-cis-decenoyl-CoA [mitochondrial matrix], ACA-CoA [mitochondrial matrix], ACAD10 [mitochondrial matrix], ACAD11 [mitochondrial inner membrane], ACADL tetramer [mitochondrial matrix], ACADM tetramer [mitochondrial matrix], ACADS tetramer [mitochondrial matrix] .... Reactions: Linoleoyl-CoA is converted to cis,cis-3,6- dodecadienoyl-CoA by three cycles of beta-oxidation; dehydrogenation of 4-cis-decenoyl-CoA to form 2-trans-4-cis-decadienoyl-CoA; Trans,cis-lauro-2,6-dienoyl-CoA is converted to cis,cis-3,6- 4-cis-decenoyl-CoA by one cycle of beta-oxidation; Reduction of 2-trans-4-cis-decadienoyl-CoA to form 3-trans-decenoyl-CoA ....
2. **EC 2.3.1 (transferases) / EC 2.3.1.16** - level 2, basis function, identifier `ec-code:2.3.1.16`, 24 nodes (17 species, 7 reactions). Species: 3-Oxodecanoyl-CoA [mitochondrial matrix], 3-Oxododecanoyl-CoA [mitochondrial matrix], 3-Oxohexanoyl-CoA [mitochondrial matrix], 3-Oxooctanoyl-CoA [mitochondrial matrix], 3-Oxopalmitoyl-CoA [mitochondrial matrix], 3-Oxotetradecanoyl-CoA [mitochondrial matrix], 3OA-CoA [mitochondrial matrix], ACAA2 tetramer [mitochondrial matrix] .... Reactions: 3-Oxotetradecanoyl-CoA+CoA-SH => Lauroyl-CoA + acetyl-CoA; 3-Oxopalmitoyl-CoA+CoA-SH => myristoyl-CoA+acetyl-CoA; 3-Oxododecanoyl-CoA+CoA-SH => Decanoyl-CoA + acetyl-CoA; 3-Oxohexanoyl-CoA+CoA-SH => Butanoyl-CoA + acetyl-CoA ....
3. **GO:0006637 / cytosol** - level 2, basis localisation, identifier `compartment_70101`, 12 nodes (10 species, 2 reactions). Species: ACBD6 [cytosol], ACBD6:LCFA-CoA [cytosol], CoA-SH [cytosol], Cytosolic ACOTs, THEM4 [cytosol], H+ [cytosol], H2O [cytosol], LCFA-CoA [cytosol], MCFA, LCFA [cytosol] .... Reactions: Cytosolic ACOTs hydrolyse MCFA-CoA, LCFA-CoA; ACBD6 binds LCFA-CoA.

## skm_pss

- Nodes: 1817 (961 species, 856 reactions). Groups: 98 (L1 12, L2 34, L3 52).
- Size distribution (nodes per group): L1: n=12, min 7, median 91, max 473, groups under 5 nodes 0; L2: n=34, min 6, median 24, max 313, groups under 5 nodes 0; L3: n=52, min 3, median 9, max 189, groups under 5 nodes 10.
- Unnamed groups: 0 of 98 (0%). Name sources: compartment 52, notes 44, sbo 1, xref 1. Labels that are a bare identifier (GO id, pathway id, EC number without a name): 1.
- Largest level-1 box: 0.26 of all nodes.
- Reactions by rule: R1 20, R10 2, R2 403, R3 274, R8 157. Species by rule: R11 444, R12 517.
- Localisation: 64% of nodes sit in a localisation box (R14); 76% of nodes sit in an innermost box whose members all share one `loc_space`.
- Dropped by the page rule (R6): 19 boxes holding 118 reactions when dropped: `biocyc:PHOTOALL-PWY` (3); `kegg.pathway:map00195` (3); `biocyc:SUCSYN-PWY` (3); `kegg.pathway:map00710` (3); `kegg.pathway:map00906` (3); `mapman:2 (Cellular respiration)` (3); `SBO:0000176 (biochemical reaction)` (4); `ec-code:1 (oxidoreductases)` (4); `mapman:5 (Lipid metabolism)` (4); `mapman:30 (Clade-specific metabolism)` (4); `mapman:17 (Protein biosynthesis)` (6); `mapman:22 (Vesicle trafficking)` (6); `ec-code:2 (transferases)` (8); `ec-code:3 (hydrolases)` (8); `mapman:6 (Nucleotide metabolism)` (8); `mapman:7 (Coenzyme metabolism)` (14); `biocyc:ETHYL-PWY` (3); `mapman:24 (Solute transport)` (14); `mapman:16 (RNA homeostasis)` (17).
- Boxes dropped as too small (R5) at level 1: 17.
- Level-1 boxes (nodes, label, identifier): 473 `Phytohormone action` [mapman:10]; 336 `RNA biosynthesis` [mapman:14]; 243 `Multi-process regulation` [mapman:27]; 194 `External stimuli response` [mapman:26]; 136 `Photosynthesis` [mapman:1]; 91 `Protein modification` [mapman:18]; 88 `Secondary metabolism` [mapman:8]; 88 `Protein homeostasis` [mapman:19]; 57 `Amino acid metabolism` [mapman:4]; 57 `KEGG pathway map00908` [kegg.pathway:map00908]; 47 `Redox homeostasis` [mapman:9]; 7 `non-covalent binding` [SBO:0000177].

Example groups:

1. **Phytohormone action** - level 1, basis function, identifier `mapman:10`, 473 nodes (258 species, 215 reactions). Species: 12,13-EOT, 12-OH-JA, 12-OH-JA-Ile, 13-HPOT, 8' OH-ABA, 9-cis-10'-apo-beta-carotenal, 9-cis-beta-carotene, AAO .... Reactions: activation_ATG13A_pa_cyt (activation); activation_ATG1A_pa_cyt (activation); activation_JAM_pa_cyt (activation); rx00002 (biochemical reaction) ....
2. **Phytohormone action / strigolactone** - level 2, basis function, identifier `mapman:10.9`, 24 nodes (14 species, 10 reactions). Species: 9-cis-10'-apo-beta-carotenal, 9-cis-beta-carotene, CCD, CL, CLA, CLAMT, D14, D14|MAX2|SCF .... Reactions: rx00314 (biochemical reaction); rx00315 (biochemical reaction); rx00316 (biochemical reaction); rx00317 (biochemical reaction) ....
3. **RNA biosynthesis / DNA-binding transcriptional regulation / cytoplasm** - level 3, basis localisation, identifier `c_cyt`, 189 nodes (97 species, 92 reactions). Species: ABI4, ABI4|RGL2, ABI5, ABI5|HY5, ADH1, AREB/ABF, BGLU30, CAMTA3 .... Reactions: rx00025 (genetic production); rx00098 (genetic production); rx00111 (genetic production); rx00126 (genetic production) ....

## skm_ja

- Nodes: 199 (117 species, 82 reactions). Groups: 25 (L1 9, L2 11, L3 5).
- Size distribution (nodes per group): L1: n=9, min 8, median 13, max 66, groups under 5 nodes 0; L2: n=11, min 5, median 9, max 55, groups under 5 nodes 0; L3: n=5, min 3, median 9, max 25, groups under 5 nodes 2.
- Unnamed groups: 0 of 25 (0%). Name sources: compartment 10, ec 1, notes 14. Labels that are a bare identifier (GO id, pathway id, EC number without a name): 0.
- Largest level-1 box: 0.33 of all nodes.
- Reactions by rule: R2 48, R3 23, R8 11. Species by rule: R11 49, R12 68.
- Localisation: 55% of nodes sit in a localisation box (R14); 78% of nodes sit in an innermost box whose members all share one `loc_space`.
- Opened as catch-all (R7): `mapman:10 (Phytohormone action)`; `mapman:10.7 (Phytohormone action.jasmonic acid)`.
- Single-child chains followed (R5): `mapman:5 (Lipid metabolism)`; `mapman:7 (Coenzyme metabolism)`; `mapman:7.3 (Coenzyme metabolism.S-adenosyl methionine (SAM) cycle)`; `mapman:27 (Multi-process regulation)`; `mapman:10.7.2 (Phytohormone action.jasmonic acid.perception and signal transduction)`; `mapman:10.7.4 (Phytohormone action.jasmonic acid.transport)`.
- Boxes dropped as too small (R5) at level 1: 13.
- Level-1 boxes (nodes, label, identifier): 66 `RNA biosynthesis` [mapman:14]; 36 `Phytohormone action.jasmonic acid.biosynthesis` [mapman:10.7.1]; 33 `Phytohormone action.jasmonic acid.perception and signal transduction.receptor activities` [mapman:10.7.2.1]; 16 `Lipid metabolism.fatty acid metabolism` [mapman:5.1]; 13 `Phytohormone action.jasmonic acid.conjugation and degradation` [mapman:10.7.3]; 10 `Coenzyme metabolism.S-adenosyl methionine (SAM) cycle.S-adenosyl homocysteine hydrolase *(SAHH)` [mapman:7.3.2]; 9 `EC 1 (oxidoreductases)` [ec-code:1]; 8 `Multi-process regulation.MAP kinase cascade signalling` [mapman:27.14]; 8 `Phytohormone action.jasmonic acid.transport.jasmonic acid transporter *(JAT)` [mapman:10.7.4.1].

Example groups:

1. **RNA biosynthesis** - level 1, basis function, identifier `mapman:14`, 66 nodes (34 species, 32 reactions). Species: CLH, COI1|JA-Ile|SCF, DELLA, DELLA|JAZ, DELLA|MYC, EIN3(like), EIN3(like)|JAZ, ERF1 .... Reactions: activation_WDbHLHMYB_ca_nuc (activation); rx00018 (non-covalent binding); rx00028 (genetic production); rx00044 (degradation) ....
2. **Phytohormone action.jasmonic acid.perception and signal transduction.receptor activities / regulatory protein *(MYC2)** - level 2, basis function, identifier `mapman:10.7.2.1.2`, 10 nodes (5 species, 5 reactions). Species: MYC3, MYC3|NPR1, MYC4, MYC4|NPR1, NPR1. Reactions: rx00444 (non-covalent binding); rx00446 (non-covalent binding); rx00877 (inactivation); rx00879 (inactivation) ....
3. **RNA biosynthesis / DNA-binding transcriptional regulation / nucleus** - level 3, basis localisation, identifier `c_nuc`, 25 nodes (15 species, 10 reactions). Species: DELLA, DELLA|JAZ, DELLA|MYC, EIN3(like), EIN3(like)|JAZ, JAM, MYC2, MYC|PYL .... Reactions: rx00018 (non-covalent binding); rx00046 (genetic production); rx00159 (non-covalent binding); rx00162 (non-covalent binding) ....

## skm_sa

- Nodes: 352 (193 species, 159 reactions). Groups: 30 (L1 8, L2 8, L3 14).
- Size distribution (nodes per group): L1: n=8, min 7, median 27, max 112, groups under 5 nodes 0; L2: n=8, min 4, median 17, max 107, groups under 5 nodes 1; L3: n=14, min 6, median 9, max 72, groups under 5 nodes 0.
- Unnamed groups: 0 of 30 (0%). Name sources: compartment 16, ec 1, notes 13. Labels that are a bare identifier (GO id, pathway id, EC number without a name): 0.
- Largest level-1 box: 0.32 of all nodes.
- Reactions by rule: R1 1, R2 92, R3 41, R8 25. Species by rule: R11 92, R12 101.
- Localisation: 71% of nodes sit in a localisation box (R14); 77% of nodes sit in an innermost box whose members all share one `loc_space`.
- Single-child chains followed (R5): `mapman:7 (Coenzyme metabolism)`; `mapman:7.16 (Coenzyme metabolism.phylloquinone biosynthesis)`; `mapman:8 (Secondary metabolism)`.
- Boxes dropped as too small (R5) at level 1: 9.
- Level-1 boxes (nodes, label, identifier): 112 `External stimuli response` [mapman:26]; 81 `Phytohormone action` [mapman:10]; 77 `RNA biosynthesis` [mapman:14]; 27 `Secondary metabolism.phenolics biosynthesis` [mapman:8.2]; 26 `Amino acid metabolism` [mapman:4]; 14 `Redox homeostasis` [mapman:9]; 8 `Coenzyme metabolism.phylloquinone biosynthesis.1,4-dihydroxy-2-naphthoyl-CoA thioesterase *(MenH/DHNAT)` [mapman:7.16.6]; 7 `EC 1 (oxidoreductases)` [ec-code:1].

Example groups:

1. **External stimuli response** - level 1, basis function, identifier `mapman:26`, 112 nodes (52 species, 60 reactions). Species: ADR1, ADR1-L1, ADR1-L2, ALD1, BAK1, BAP1, BAP2, BIK1 .... Reactions: activation_EDS1_pa_cyt (activation); rx00054 (non-covalent binding); rx00087 (non-covalent binding); rx00088 (biochemical reaction) ....
2. **Phytohormone action / salicylic acid** - level 2, basis function, identifier `mapman:10.8`, 30 nodes (17 species, 13 reactions). Species: CUL3A, EDS1-WRKY18, EDS5, EIN3|NPR1, ICS, IsoChor, MYC3, MYC3|NPR1 .... Reactions: rx00069 (biochemical reaction); rx00114 (activation); rx00422 (degradation); rx00423 (degradation) ....
3. **External stimuli response / pathogen / cytoplasm** - level 3, basis localisation, identifier `c_cyt`, 72 nodes (33 species, 39 reactions). Species: ADR1, ADR1-L1, ADR1-L2, ALD1, BAK1, BAP1, BAP2, BIK1 .... Reactions: activation_EDS1_pa_cyt (activation); rx00514 (biochemical reaction); rx00524 (inactivation); rx00758 (genetic production) ....

## Sensitivity to the parameters

Group counts L1/L2/L3 and, in brackets, the largest level-1 share. Each row changes one parameter; the delivered files use the first row. These runs were written to a scratch folder, not to `out/`.

| variant | brenda_merged | e_coli_core | R-HSA-77289 | skm_pss | skm_ja | skm_sa |
|---|---|---|---|---|---|---|
| default | 7/24/0 (0.27) | 8/13/0 (0.23) | 6/6/0 (0.32) | 12/34/52 (0.26) | 9/11/5 (0.33) | 8/8/14 (0.32) |
| M_MIN=2 | 7/24/0 (0.27) | 9/14/2 (0.23) | 8/8/0 (0.30) | 12/39/72 (0.26) | 9/15/15 (0.33) | 10/17/20 (0.30) |
| M_MIN=5 | 6/19/0 (0.29) | 8/2/0 (0.23) | 4/5/0 (0.36) | 12/27/35 (0.26) | 5/12/5 (0.48) | 6/10/9 (0.32) |
| PAGE=9 | 7/24/0 (0.27) | 8/13/0 (0.23) | 6/6/0 (0.32) | 9/29/48 (0.30) | 9/11/5 (0.33) | 8/8/14 (0.32) |
| PAGE=16 | 7/24/0 (0.27) | 8/13/0 (0.23) | 6/6/0 (0.32) | 16/39/56 (0.22) | 9/11/5 (0.33) | 8/8/14 (0.32) |
| S_MAX=0.33 | 7/24/0 (0.27) | 8/13/0 (0.23) | 6/6/0 (0.32) | 12/34/52 (0.26) | 9/11/5 (0.33) | 8/8/14 (0.32) |

## Known failure cases

1. **EC boxes are not pathways.** In `brenda_merged`, `e_coli_core` and `R-HSA-77289` the only usable identifiers on reactions are EC numbers (plus flat GO terms in the Reactome file), so the boxes are chemistry classes (oxidoreductases, transferases, ...). A pathway is spread over several boxes and each box mixes steps of unrelated pathways. Against pathway-style reference labels this grouping should be expected to score poorly on these three models. The consequence is visible in the clones: 44 of 178 species in `brenda_merged` and 25 of 72 in `e_coli_core` touch three or more level-1 boxes.
2. **Many labels are bare identifiers.** The bundle has no EC names and no GO names. Below the seven EC classes a label is the number alone (`EC 2.7`), a GO box is `GO:0006637`, a KEGG box is `KEGG pathway map00908`. Counts of such labels: brenda_merged 24 of 31, e_coli_core 11 of 21, R-HSA-77289 4 of 12, skm_pss 1 of 98, skm_ja 0 of 25, skm_sa 0 of 30. They are backed by an identifier but a reader has to look the identifier up. In the SKM models the EC names come from the MapMan notes, so this does not occur there.
3. **No depth can be chosen for GO and pathway cross-references.** There is no GO graph in `common/`, so a GO box is one exact term and biological-process and molecular-function terms cannot be told apart. In `R-HSA-77289` most GO terms repeat the EC number of the same reaction; GO is used only for the 4 reactions EC does not place (box `GO:0006637`).
4. **The page rule moves reactions into boxes that do not describe them (`skm_pss`).** 19 boxes with 118 reactions (14% of the reactions) were dropped to keep 12 level-1 boxes, among them `Solute transport`, `Coenzyme metabolism` and `RNA homeostasis`. Their reactions are then placed by propagation. Example: the three SAMS species (own bin: Coenzyme metabolism, SAM cycle) end up in `Amino acid metabolism` (1) and `Protein homeostasis` (2).
5. **Propagation puts unannotated nodes under a label that was not assigned to them.** Share of reactions placed by neighbours (R8, R10): brenda_merged 0%, e_coli_core 1%, R-HSA-77289 4%, skm_pss 19%, skm_ja 13%, skm_sa 16%. All metabolites and complexes (they carry no usable term) are placed by their reactions (R12). The label of a box is therefore backed by an identifier only for its annotated members.
6. **Ambiguous proteins are decided by tie-breaks.** In the SKM models most proteins carry several MapMan bins (reactions whose bins tie between several top bins after R4b: `skm_ja` 49 of 76, `skm_sa` 33 of 136, `skm_pss` 263 of 743). R4b-d decide by sources, neighbours and box size; a different but equally defensible reading exists for many of them. Example: in `skm_ja` the JA regulator MYC2 (bins: jasmonic acid signalling, and bHLH transcription factor) is placed in `RNA biosynthesis`, while MYC3 and MYC4 with the same bins are in the jasmonic acid box.
7. **Family bins act as near catch-alls.** `RNA biosynthesis` (transcription factors and the gene-expression reactions they regulate) holds 66 of 199 nodes in `skm_ja`, including targets from unrelated pathways and, through R12, nodes such as the complex `COI1|JA-Ile|SCF`. It stays below `S_MAX`, so it is not opened. In `skm_pss` its sub-box `DNA-binding transcriptional regulation` holds 313 of its 336 nodes.
8. **Uneven depth and long labels at level 1.** After a single-child chain (R5) a small box is named by a very specific term, next to broad boxes: in `skm_ja`, `Coenzyme metabolism.S-adenosyl methionine (SAM) cycle.S-adenosyl homocysteine hydrolase *(SAHH)` (10 nodes) stands beside `RNA biosynthesis` (66). The label is exact for the annotated members but it is a protein-level box.
9. **Only two function levels.** Level 2 is one tree level below level 1. In `skm_pss` the box `Phytohormone action / jasmonic acid` (111 nodes) is not divided into biosynthesis and signalling, although the bins would allow it; in `skm_ja`, where R7 opens two levels, it is.
10. **Dropping is greedy.** Boxes that are too small are dropped one at a time in a fixed order (R5) and are never re-admitted, so a box dropped early might have reached 3 reactions after a later drop. The order (ontology, depth, potential size, size) limits this but does not exclude it.
11. **Function boxes span compartments.** A localisation part is made only for places with at least 3 nodes and only if two places qualify, to avoid tiny boxes. Nodes of small places stay loose in the function box: skm_pss 24%, skm_ja 22%, skm_sa 23% of nodes are in an innermost box that mixes `loc_space` values (0% in the three metabolic models). Level-1 boxes that have function sub-boxes are never divided by place for their loose nodes.
12. **Leftovers.** `skm_pss`: two binding reactions and five species (GPAphid2, RANGAP, Rx and their complexes) are a component no term reaches; they form a 7-node level-1 box named by their common SBO term `non-covalent binding`. `e_coli_core`: the biomass reaction is alone in its SBO class and is placed by R10 with most of its participants (an EC box that does not describe it). `R-HSA-77289`: 2 reactions placed by R10.
13. **Two reading rules are tailored to conventions, not to a standard.** (a) A bin whose names are `EC_<n> ...` is read as the EC tree; this is how MapMan writes its bin 50, and it is what keeps `Enzyme classification` from absorbing every enzyme, but another bin vocabulary would not trigger it. (b) Pathway cross-references are recognised by identifier syntax (KEGG map numbers, `PWY` frames); a pathway identifier in another syntax is not seen. The seven EC class names are embedded in the code.
14. **How the rules were arrived at.** `M_MIN`, `PAGE` and `S_MAX` were fixed before the first run. The tie-break R4c-d, the EC-in-bin reading, the lone-sub-box condition of R13, the two-place condition of R14 and the unit of the clone rule were revised after looking at the output on the six panel models (judged by eye on known proteins, without reference labels). They are general rules, but they have not been tried on any model outside this panel.
