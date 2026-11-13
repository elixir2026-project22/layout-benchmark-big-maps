# NOTES - strategy `location_type`

Cell localisation on top (level 1), reaction type inside each compartment (level 2), one more type
level only inside a box of more than 50 nodes (level 3). The rules are in `RULES.md`.

The overview and the per-model sections were generated from the delivered `out/<model>/groups.tsv` and
`groups_meta.tsv` and the feature tables of the bundle. The numbers quoted in the failure list were
read from the same files. "Size" of a group counts all nodes under it,
nested members included.

## Overview (delivered variant: `out/`, single-compartment level skipped, rule R3)

| model | L1 | L2 | L3 | unnamed | largest L1 share | groups < 5 nodes | innermost boxes > 50 nodes | clone copies | validate | two runs identical | groups.tsv sha256 (first 16) | groups_meta.tsv sha256 (first 16) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| brenda_merged | 7 | 17 | 0 | 0 | 0.288 | 0 | 0 | 211 | valid | yes | `cd1db0968cad6e54` | `8741ee1533379fc8` |
| e_coli_core | 3 | 6 | 0 | 0 | 0.611 | 0 | 0 | 41 | valid | yes | `262794a099167bd6` | `c13ef2c8e0173d0d` |
| R-HSA-77289 | 4 | 5 | 2 | 0 | 0.865 | 2 | 0 | 9 | valid | yes | `c41ab0e6d10ba1fc` | `18d79198b177e549` |
| skm_pss | 17 | 18 | 4 | 0 | 0.733 | 5 | 7 | 124 | valid | yes | `0a8263a6ea14df1f` | `5911249932a63956` |
| skm_ja | 9 | 8 | 0 | 0 | 0.487 | 4 | 0 | 8 | valid | yes | `98983a3aaf21bac4` | `c801c11d4e6903fd` |
| skm_sa | 10 | 12 | 0 | 0 | 0.551 | 3 | 2 | 21 | valid | yes | `e6e5cdfb99ca4a87` | `f94a3d4de5212db4` |

- `python common/contract.py validate` printed `valid` for all six models.
- Determinism: `strategy.py` was run twice per model into separate folders; the sha256 of `groups.tsv`
  and of `groups_meta.tsv` is the same in both runs for all six models (hashes in the table). A third
  run used the PowerShell loop of the docstring and gave the same hashes for all six models. That run was
  started as a PowerShell subprocess with the Python interpreter's folder added to PATH: the sandbox's
  own PowerShell tool has no `python` on PATH, so a first attempt there wrote no files (not a
  difference in output).
- The `*_grouped.xml` files were written with `python common/contract.py write-sbml`; each holds as many
  SBML groups as `groups_meta.tsv` has rows.
- libsbml reports one error when reading `R-HSA-77289_grouped.xml` (a `<layout>` without `layout:id`).
  The same error is reported for the input `panel/models/R-HSA-77289.xml`; it is carried over, not introduced.

## Per model

### brenda_merged

396 nodes (218 reactions, 178 species), 24 groups.

| level | groups | size min / median / max (nodes, nested members included) |
|---|---|---|
| L1 | 7 | 8 / 47 / 114 |
| L2 | 17 | 6 / 12 / 41 |
| L3 | 0 | – |

- Size distribution over all groups (nodes per group): 1: 0, 2–4: 0, 5–10: 8, 11–25: 8, 26–50: 5, 51–100: 1, >100: 2.
- Unnamed groups: 0 of 24 (0%). Name sources: compartment 1, ec 23.
- Largest level-1 group: EC 1 oxidoreductases, 28.8% of the nodes.
- Groups under 5 nodes: 0 (type boxes among them: 0); singletons: 0.
- Innermost boxes over 50 nodes: 0; the largest holds 47 nodes (`ec_4`).
- Nodes that sit directly in a group that has sub-boxes: 5 (3 reactions, 2 species).
- Clones: 211 copies of 62 species.
- Participant edges whose species is homed in another box than the reaction: 472 of 914 (52%); 78 (9%) once the clones are drawn.
- Nodes per rule: R5 65, R6 149, R8 4, R9 178.

Three example groups (chosen by a fixed recipe: the largest boundary group, or the smallest localisation group where the model has no boundary group; the first-type-level box closest to the median size; the largest box of the deepest level):

- **default** — `loc_c`, level 1, named from `compartment` `c`; 8 nodes = 4 reactions + 4 species (the smallest localisation group).
  - reactions (4 of 4): R_rx_NODE_4393 (NODE_4393); R_rx_NODE_4439 (NODE_4439); R_rx_NODE_4449 (NODE_4449); R_rx_NODE_4455 (NODE_4455)
  - species (4 of 4): oxidized electron donor; reduced electron acceptor; 5-ureido-4-imidazole carboxylate; 4-aminoimidazole
- **EC 3 hydrolases** — `ec_3`, level 1, named from `ec` `ec-code:3`; 55 nodes = 35 reactions + 20 species (the first-type-level box closest to the median size).
  - reactions (5 of 35): R_rx_0b67f0c415 (BR22769); R_rx_0fb949fb5f (BR22743); R_rx_28809f586e (BR16044); R_rx_288676c879 (BR15958); R_rx_359808fee2 (BR22683)
  - species (5 of 20): dGTP; 5-amino-4-imidazole carboxylate; IMP; H2O; Xanthosine
- **EC 2 transferases / EC 2.7** — `ec_2/ec_2.7`, level 2, named from `ec` `ec-code:2.7`; 41 nodes = 26 reactions + 15 species (the largest box of the deepest level).
  - reactions (5 of 26): R_rx_189c95a8a7 (BR47825); R_rx_1a1985edf3 (BR48070); R_rx_2bff8e5db7 (BS353783); R_rx_35de376d34 (BR47850); R_rx_3db9ac4c4f (BR47726)
  - species (5 of 15): phosphoenolpyruvate; GDP; dAMP; GMP; ADP

### e_coli_core

167 nodes (95 reactions, 72 species), 9 groups.

| level | groups | size min / median / max (nodes, nested members included) |
|---|---|---|
| L1 | 3 | 25 / 40 / 102 |
| L2 | 6 | 5 / 11.5 / 32 |
| L3 | 0 | – |

- Size distribution over all groups (nodes per group): 1: 0, 2–4: 0, 5–10: 3, 11–25: 2, 26–50: 3, 51–100: 0, >100: 1.
- Unnamed groups: 0 of 9 (0%). Name sources: compartment 3, ec 6.
- Largest level-1 group: cytosol, 61.1% of the nodes.
- Groups under 5 nodes: 0 (type boxes among them: 0); singletons: 0.
- Innermost boxes over 50 nodes: 0; the largest holds 40 nodes (`loc_e`).
- Nodes that sit directly in a group that has sub-boxes: 8 (3 reactions, 5 species).
- Clones: 41 copies of 20 species.
- Participant edges whose species is homed in another box than the reaction: 198 of 360 (55%); 137 (38%) once the clones are drawn.
- Nodes per rule: R1 40, R2 25, R5 47, R8 3, R9 51, R10 1.

Three example groups (chosen by a fixed recipe: the largest boundary group, or the smallest localisation group where the model has no boundary group; the first-type-level box closest to the median size; the largest box of the deepest level):

- **cytosol boundary** — `loc_c@boundary`, level 1, named from `compartment` `c`; 25 nodes = 25 reactions + 0 species (the largest boundary group).
  - reactions (5 of 25): R_ACALDt (Acetaldehyde reversible transport); R_ACt2r (Acetate reversible transport via proton symport); R_AKGt2r (2 oxoglutarate reversible transport via symport); R_ATPS4r (ATP synthase (four protons for one ATP)); R_CO2t (CO2 transporter via diffusion)
- **cytosol / EC 4 lyases** — `loc_c/ec_4`, level 2, named from `ec` `ec-code:4`; 16 nodes = 8 reactions + 8 species (the first-type-level box closest to the median size).
  - reactions (5 of 8): R_ACONTa (Aconitase (half-reaction A, Citrate hydro-lyase)); R_ACONTb (Aconitase (half-reaction B, Isocitrate hydro-lyase)); R_ENO (Enolase); R_FBA (Fructose-bisphosphate aldolase); R_FUM (Fumarase)
  - species (5 of 8): D-Glycerate 2-phosphate; Cis-Aconitate; Dihydroxyacetone phosphate; Glyoxylate; H2O H2O
- **cytosol / EC 1 oxidoreductases** — `loc_c/ec_1`, level 2, named from `ec` `ec-code:1`; 32 nodes = 15 reactions + 17 species (the largest box of the deepest level).
  - reactions (5 of 15): R_ACALD (Acetaldehyde dehydrogenase (acetylating)); R_AKGDH (2-Oxogluterate dehydrogenase); R_ALCD2x (Alcohol dehydrogenase (ethanol)); R_G6PDH2r (Glucose 6-phosphate dehydrogenase); R_GAPD (Glyceraldehyde-3-phosphate dehydrogenase)
  - species (5 of 17): 3-Phospho-D-glyceroyl phosphate; 6-phospho-D-glucono-1,5-lactone; Acetaldehyde; 2-Oxoglutarate; CO2 CO2

### R-HSA-77289

148 nodes (47 reactions, 101 species), 11 groups.

| level | groups | size min / median / max (nodes, nested members included) |
|---|---|---|
| L1 | 4 | 4 / 8 / 128 |
| L2 | 5 | 10 / 17 / 56 |
| L3 | 2 | 17 / 26.5 / 36 |

- Size distribution over all groups (nodes per group): 1: 0, 2–4: 2, 5–10: 1, 11–25: 5, 26–50: 1, 51–100: 1, >100: 1.
- Unnamed groups: 0 of 11 (0%). Name sources: compartment 4, ec 7.
- Largest level-1 group: mitochondrial matrix, 86.5% of the nodes.
- Groups under 5 nodes: 2 (type boxes among them: 0); singletons: 0.
- Innermost boxes over 50 nodes: 0; the largest holds 36 nodes (`loc_compartment_5460/ec_1/ec_1.3`).
- Nodes that sit directly in a group that has sub-boxes: 14 (4 reactions, 10 species).
- Clones: 9 copies of 5 species.
- Participant edges whose species is homed in another box than the reaction: 78 of 256 (30%); 60 (23%) once the clones are drawn.
- Nodes per rule: R1 16, R2 4, R5 22, R6 18, R8 3, R9 83, R10 2.

Three example groups (chosen by a fixed recipe: the largest boundary group, or the smallest localisation group where the model has no boundary group; the first-type-level box closest to the median size; the largest box of the deepest level):

- **mitochondrial inner membrane** — `loc_compartment_5460@boundary`, level 1, named from `compartment` `compartment_70616`; 4 nodes = 1 reactions + 3 species (the largest boundary group).
  - reactions (1 of 1): reaction_8933547 (MCAT transfers Mal from Mal-CoA to NDUFAB1)
  - species (3 of 3): NDUFAB1 [mitochondrial inner membrane]; ACAD11 [mitochondrial inner membrane]; Mal-PPANT-S112-NDUFAB1 [mitochondrial inner membrane]
- **mitochondrial matrix / EC 6 ligases** — `loc_compartment_5460/ec_6`, level 2, named from `ec` `ec-code:6`; 17 nodes = 3 reactions + 14 species (the first-type-level box closest to the median size).
  - reactions (3 of 3): reaction_5696004 (ACSF2 ligates CoA-SH to MCFA); reaction_71031 (propionyl-CoA + CO2 + ATP <=> D-methylmalonyl-CoA + ADP + orthophosphate); reaction_8875013 (ACSM3,ACSM6 ligate CoA to BUT)
  - species (5 of 14): Pi [mitochondrial matrix]; ADP [mitochondrial matrix]; ATP [mitochondrial matrix]; AMP [mitochondrial matrix]; PPi [mitochondrial matrix]
- **mitochondrial matrix / EC 1 oxidoreductases / EC 1.3** — `loc_compartment_5460/ec_1/ec_1.3`, level 3, named from `ec` `ec-code:1.3`; 36 nodes = 11 reactions + 25 species (the largest box of the deepest level).
  - reactions (5 of 11): reaction_109341 (dehydrogenation of 4-cis-decenoyl-CoA to form 2-trans-4-cis-decadienoyl-CoA); reaction_5695980 (ACAD10 dehydrogenates S-2MPDA-CoA); reaction_5695989 (ACAD11 dehydrogenates BH-CoA); reaction_77263 (lauroyl-CoA+FAD => 2-trans-Dodecenoyl-CoA+FADH2); reaction_77274 (myristoyl-CoA+FAD => trans-tetradec-2-enoyl-CoA+FADH2)
  - species (5 of 25): 2-trans-4-cis-decadienoyl-CoA [mitochondrial matrix]; H+ [mitochondrial matrix]; NADP+ [mitochondrial matrix]; FAD [mitochondrial matrix]; NADPH [mitochondrial matrix]

### skm_pss

1817 nodes (856 reactions, 961 species), 39 groups.

| level | groups | size min / median / max (nodes, nested members included) |
|---|---|---|
| L1 | 17 | 1 / 10 / 1331 |
| L2 | 18 | 7 / 39 / 395 |
| L3 | 4 | 14 / 37.5 / 240 |

- Size distribution over all groups (nodes per group): 1: 1, 2–4: 4, 5–10: 6, 11–25: 9, 26–50: 6, 51–100: 6, >100: 7.
- Unnamed groups: 0 of 39 (0%). Name sources: compartment 17, sbo 22.
- Largest level-1 group: cytoplasm, 73.3% of the nodes.
- Groups under 5 nodes: 5 (type boxes among them: 0); singletons: 1.
- Innermost boxes over 50 nodes: 7; the largest holds 395 nodes (`loc_c_cyt/sbo_0000177`).
- Nodes that sit directly in a group that has sub-boxes: 67 (6 reactions, 61 species).
- Clones: 124 copies of 87 species.
- Participant edges whose species is homed in another box than the reaction: 698 of 2257 (31%); 558 (25%) once the clones are drawn.
- Nodes per rule: R1 68, R2 46, R4 589, R7 190, R8 6, R9 859, R10 59.

Three example groups (chosen by a fixed recipe: the largest boundary group, or the smallest localisation group where the model has no boundary group; the first-type-level box closest to the median size; the largest box of the deepest level):

- **nucleus boundary** — `loc_c_nuc@boundary`, level 1, named from `compartment` `c_nuc`; 64 nodes = 64 reactions + 0 species (the largest boundary group).
  - reactions (5 of 64): rx00016 (EBF → EBF|SCF); rx00043 (COI1 + JA-Ile → COI1|JA-Ile|SCF); rx00090 (EDS1|PAD4 + MPK3,6 → EDS1|MPK3|PAD4); rx00097 (NPR1 → NPR1); rx00122 (WRKY53 + ESP → ESP|WRKY53)
- **chloroplast / biochemical reaction** — `loc_c_chl/sbo_0000176`, level 2, named from `sbo` `SBO:0000176`; 35 nodes = 11 reactions + 24 species (the first-type-level box closest to the median size).
  - reactions (5 of 11): rx00014 (ROS + ALA → MDA); rx00026 (ALA → 13-HPOT); rx00029 (13-HPOT → 12,13-EOT); rx00030 (12,13-EOT → OPDA); rx00069 (Chor → IsoChor)
  - species (5 of 24): 12,13-EOT; 13-HPOT; ADT; ALA; AOC
- **cytoplasm / genetic production / with stimulator** — `loc_c_cyt/sbo_0000589/role_0000459`, level 3, named from `sbo` `SBO:0000459`; 240 nodes = 144 reactions + 96 species (the largest box of the deepest level).
  - reactions (5 of 144): rx00025 (∅ → PDF1.2); rx00028 (∅ → LOX); rx00049 (∅ → CLH); rx00050 (∅ → JR1); rx00051 (∅ → PR13)
  - species (5 of 96): ABI4|RGL2; ABI4; ABI5; ADH1; ADR1-L1

### skm_ja

199 nodes (82 reactions, 117 species), 17 groups.

| level | groups | size min / median / max (nodes, nested members included) |
|---|---|---|
| L1 | 9 | 1 / 9 / 97 |
| L2 | 8 | 5 / 13.5 / 23 |
| L3 | 0 | – |

- Size distribution over all groups (nodes per group): 1: 1, 2–4: 3, 5–10: 3, 11–25: 7, 26–50: 2, 51–100: 1, >100: 0.
- Unnamed groups: 0 of 17 (0%). Name sources: compartment 9, sbo 8.
- Largest level-1 group: cytoplasm, 48.7% of the nodes.
- Groups under 5 nodes: 4 (type boxes among them: 0); singletons: 1.
- Innermost boxes over 50 nodes: 0; the largest holds 30 nodes (`loc_c_per`).
- Nodes that sit directly in a group that has sub-boxes: 20 (2 reactions, 18 species).
- Clones: 8 copies of 5 species.
- Participant edges whose species is homed in another box than the reaction: 77 of 218 (35%); 67 (31%) once the clones are drawn.
- Nodes per rule: R1 52, R2 16, R4 47, R8 2, R9 65, R10 17.

Three example groups (chosen by a fixed recipe: the largest boundary group, or the smallest localisation group where the model has no boundary group; the first-type-level box closest to the median size; the largest box of the deepest level):

- **nucleus boundary** — `loc_c_nuc@boundary`, level 1, named from `compartment` `c_nuc`; 9 nodes = 9 reactions + 0 species (the largest boundary group).
  - reactions (5 of 9): rx00043 (COI1 + JA-Ile → COI1|JA-Ile|SCF); rx00128 (WRKY57 + JAZ → JAZ|WRKY57); rx00232 (RBC + COI1 → COI1|RBC); rx00233 (RBC + JAZ → JAZ|RBC); rx00459 (MYC2 + PYL → MYC|PYL)
- **cytoplasm / non-covalent binding** — `loc_c_cyt/sbo_0000177`, level 2, named from `sbo` `SBO:0000177`; 13 nodes = 4 reactions + 9 species (the first-type-level box closest to the median size).
  - reactions (4 of 4): rx00444 (NPR1 + MYC3 → MYC3|NPR1); rx00446 (NPR1 + MYC4 → MYC4|NPR1); rx00634 (VPg + LOX → LOX2|VPg); rx00679 (RCA + VPg → RCA|VPg)
  - species (5 of 9): LOX2|VPg; MYC3|NPR1; MYC3; MYC4|NPR1; MYC4
- **cytoplasm / biochemical reaction** — `loc_c_cyt/sbo_0000176`, level 2, named from `sbo` `SBO:0000176`; 23 nodes = 7 reactions + 16 species (the largest box of the deepest level).
  - reactions (5 of 7): rx00020 (Thr → Ile); rx00040 (JA → MeJA); rx00041 (Ile + JA → JA-Ile); rx00042 (JA-Ile → 12-OH-JA-Ile); rx00594 (SAH → L-homo-cys)
  - species (5 of 16): 12-OH-JA-Ile; 12-OH-JA; CYP94; IAR; Ile

### skm_sa

352 nodes (159 reactions, 193 species), 22 groups.

| level | groups | size min / median / max (nodes, nested members included) |
|---|---|---|
| L1 | 10 | 1 / 10 / 194 |
| L2 | 12 | 5 / 12.5 / 84 |
| L3 | 0 | – |

- Size distribution over all groups (nodes per group): 1: 3, 2–4: 0, 5–10: 7, 11–25: 4, 26–50: 4, 51–100: 3, >100: 1.
- Unnamed groups: 0 of 22 (0%). Name sources: compartment 10, sbo 12.
- Largest level-1 group: cytoplasm, 55.1% of the nodes.
- Groups under 5 nodes: 3 (type boxes among them: 0); singletons: 3.
- Innermost boxes over 50 nodes: 2; the largest holds 84 nodes (`loc_c_cyt/sbo_0000589`).
- Nodes that sit directly in a group that has sub-boxes: 24 (1 reactions, 23 species).
- Clones: 21 copies of 13 species.
- Participant edges whose species is homed in another box than the reaction: 174 of 404 (43%); 149 (37%) once the clones are drawn.
- Nodes per rule: R1 32, R2 14, R4 134, R8 1, R9 148, R10 23.

Three example groups (chosen by a fixed recipe: the largest boundary group, or the smallest localisation group where the model has no boundary group; the first-type-level box closest to the median size; the largest box of the deepest level):

- **nucleus boundary** — `loc_c_nuc@boundary`, level 1, named from `compartment` `c_nuc`; 13 nodes = 13 reactions + 0 species (the largest boundary group).
  - reactions (5 of 13): rx00090 (EDS1|PAD4 + MPK3,6 → EDS1|MPK3|PAD4); rx00097 (NPR1 → NPR1); rx00500 (NPR1 + WRKY6,18 → NPR1|WRKY18); rx00501 (NPR1 + CDK → CDK|NPR1); rx00502 (WRKY6,18 + CDK → CDK|WRKY6,18)
- **chloroplast / biochemical reaction** — `loc_c_chl/sbo_0000176`, level 2, named from `sbo` `SBO:0000176`; 13 nodes = 4 reactions + 9 species (the first-type-level box closest to the median size).
  - reactions (4 of 4): rx00069 (Chor → IsoChor); rx00073 (Chor → Prep); rx00076 (Phe → CA); rx00541 (L-arogenate → Phe)
  - species (5 of 9): ADT; CA; CM; Chor; IsoChor
- **cytoplasm / genetic production** — `loc_c_cyt/sbo_0000589`, level 2, named from `sbo` `SBO:0000589`; 84 nodes = 48 reactions + 36 species (the largest box of the deepest level).
  - reactions (5 of 48): rx00098 (∅ → WRKY70); rx00101 (∅ → PR1); rx00102 (∅ → PR2); rx00103 (∅ → PR5); rx00148 (∅ → EDS1)
  - species (5 of 36): ADR1-L1; ADR1-L2; ADR1; ALD1; BAK1

## Variant `--keep-single-compartment` (R3 switched off)

Folder `out_keep_single/`. For the five models with more than one compartment the files are
byte-identical to `out/` (checked by hash: e_coli_core, R-HSA-77289, skm_pss, skm_ja, skm_sa), so only `brenda_merged` is
included there. validate: `valid`; two runs identical: yes
(groups.tsv `10994ef6a562a16e`, groups_meta.tsv `2d05df0712a5ccd7`).

### brenda_merged, variant keep_single_compartment

396 nodes (218 reactions, 178 species), 24 groups.

| level | groups | size min / median / max (nodes, nested members included) |
|---|---|---|
| L1 | 1 | 396 / 396 / 396 |
| L2 | 6 | 20 / 51 / 114 |
| L3 | 17 | 6 / 12 / 41 |

- Size distribution over all groups (nodes per group): 1: 0, 2–4: 0, 5–10: 7, 11–25: 8, 26–50: 5, 51–100: 1, >100: 3.
- Unnamed groups: 0 of 24 (0%). Name sources: compartment 1, ec 23.
- Largest level-1 group: default, 100.0% of the nodes.
- Groups under 5 nodes: 0 (type boxes among them: 0); singletons: 0.
- Innermost boxes over 50 nodes: 0; the largest holds 47 nodes (`loc_c/ec_4`).
- Nodes that sit directly in a group that has sub-boxes: 13 (7 reactions, 6 species).
- Clones: 208 copies of 62 species.
- Participant edges whose species is homed in another box than the reaction: 472 of 914 (52%); 81 (9%) once the clones are drawn.
- Nodes per rule: R5 65, R6 149, R8 4, R9 178.

Three example groups (chosen by a fixed recipe: the largest boundary group, or the smallest localisation group where the model has no boundary group; the first-type-level box closest to the median size; the largest box of the deepest level):

- **default** — `loc_c`, level 1, named from `compartment` `c`; 396 nodes = 218 reactions + 178 species (the smallest localisation group).
  - reactions (5 of 218): R_rx_028c1c2a2f (BR47671); R_rx_058ab7bb35 (BR43957); R_rx_08b3253feb (BR50094); R_rx_08c4822df6 (BS287796); R_rx_0aac87f3c2 (BR53166)
  - species (5 of 178): dGTP; purine; 10-formyltetrahydrofolate; (3S)-citryl-CoA; 5-amino-4-imidazole carboxylate
- **default / EC 3 hydrolases** — `loc_c/ec_3`, level 2, named from `ec` `ec-code:3`; 55 nodes = 35 reactions + 20 species (the first-type-level box closest to the median size).
  - reactions (5 of 35): R_rx_0b67f0c415 (BR22769); R_rx_0fb949fb5f (BR22743); R_rx_28809f586e (BR16044); R_rx_288676c879 (BR15958); R_rx_359808fee2 (BR22683)
  - species (5 of 20): dGTP; 5-amino-4-imidazole carboxylate; IMP; H2O; Xanthosine
- **default / EC 2 transferases / EC 2.7** — `loc_c/ec_2/ec_2.7`, level 3, named from `ec` `ec-code:2.7`; 41 nodes = 26 reactions + 15 species (the largest box of the deepest level).
  - reactions (5 of 26): R_rx_189c95a8a7 (BR47825); R_rx_1a1985edf3 (BR48070); R_rx_2bff8e5db7 (BS353783); R_rx_35de376d34 (BR47850); R_rx_3db9ac4c4f (BR47726)
  - species (5 of 15): phosphoenolpyruvate; GDP; dAMP; GMP; ADP

The same 23 EC boxes, one level deeper, under a single level-1 box named "default" (the name the SBML
gives its one compartment) that holds 100% of the model. Recommended is the delivered variant in
`out/`: a level-1 box that holds everything tells the reader nothing, and the contract lists "one
compartment only" among the cases that need a fallback.

## Known failure cases

1. **One compartment holds most of the model.** Level 1 is given by the model's localisation, so the
   strategy cannot avoid a catch-all where the biology sits in one compartment. Largest level-1 share:
   brenda_merged 28.8% (EC 1 oxidoreductases), e_coli_core 61.1% (cytosol), R-HSA-77289 86.5% (mitochondrial matrix), skm_pss 73.3% (cytoplasm), skm_ja 48.7% (cytoplasm), skm_sa 55.1% (cytoplasm).
   In four of the six models one box holds more than half of the nodes.
2. **Boxes that location and type together cannot divide.** After level 3, these innermost boxes still
   hold more than 50 nodes (MAX_BOX):
   - brenda_merged: none
   - e_coli_core: none
   - R-HSA-77289: none
   - skm_pss: cytoplasm / non-covalent binding (395); cytoplasm / biochemical reaction (308); cytoplasm / activation (253); cytoplasm / genetic production / with stimulator (240); nucleus / non-covalent binding (52); nucleus / genetic production / with stimulator (52); nucleus / activation (51)
   - skm_ja: none
   - skm_sa: cytoplasm / genetic production (84); cytoplasm / biochemical reaction (61)
   In the Stress Knowledge Map models the reactions carry one SBO term each and almost no EC number
   (4 of 856 reactions in skm_pss, 0 of 82 in skm_ja, 3 of 159 in skm_sa), and within one SBO term the
   modifier role and the entity type of the products are nearly constant, so no type attribute in the
   SBML separates these boxes further. Dividing them needs function or topology, which belong to other tracks.
3. **Level 3 by modifier role is lopsided.** In skm_pss, cytoplasm / genetic production splits into
   "with stimulator" (240 nodes) and "with inhibitor" (23 nodes); nucleus / genetic production into 52 and 14.
   The split is meaningful (induced against repressed expression) but relieves the large box very little.
4. **Tiny level-1 groups.** A compartment or boundary with few nodes is still its own level-1 group,
   because localisation is given (contract rule 5). Groups under 5 nodes:
   - brenda_merged: none
   - e_coli_core: none
   - R-HSA-77289: endoplasmic reticulum lumen (4); mitochondrial inner membrane (4)
   - skm_pss: golgi apparatus (1); mitochondrion (3); mitochondrion boundary (3); vacuole (4); vacuole boundary (3)
   - skm_ja: chloroplast boundary (2); endoplasmic reticulum (4); endoplasmic reticulum boundary (1); peroxisome boundary (4)
   - skm_sa: endoplasmic reticulum (1); vacuole (1); vacuole boundary (1)
5. **Type boxes cut pathways.** A pathway runs through several reaction types, so its metabolites and
   proteins are produced in one box and used in another. Share of participant edges that leave the
   reaction's box (before / after the clones are drawn): brenda_merged 52% / 9%, e_coli_core 55% / 38%, R-HSA-77289 30% / 23%, skm_pss 31% / 25%, skm_ja 35% / 31%, skm_sa 43% / 37%.
   These counts include the edges of transport reactions, which always cross (the reaction sits in the
   boundary group, its participants inside the two spaces).
6. **Boundary groups hold reactions only.** Species are located inside a space, never on its boundary
   (except the three membrane species of R-HSA-77289), so a boundary box shows transport steps without
   their participants; e.g. "cytosol boundary" of e_coli_core is 25 reactions and no species.
7. **Nodes outside the type boxes.** Reactions of a rare type and species that take part only in
   transport sit directly in the compartment box, next to its type boxes: brenda_merged 5, e_coli_core 8, R-HSA-77289 14, skm_pss 67, skm_ja 20, skm_sa 24 nodes.
   In skm_pss 59 of these are species with no reaction inside their own level-1 group (R10); 38 of the 59
   take part in boundary (transport) reactions only, the other 21 also in reactions of another compartment.
8. **Labels that say little.** EC sub-classes are labelled by number only ("EC 1.3"), because the
   strategy carries the names of the seven EC classes and nothing below them. "biochemical reaction" is
   the SBO term itself where the model uses it, but it also receives rolled-up rare terms (R4), so it is
   the least specific box of a compartment. "default" is the SBML name of the one compartment of
   brenda_merged. Level-3 labels are the SBO role name with "with" in front ("with stimulator").
9. **Reactions with EC numbers of several classes** are assigned to the class most of their EC numbers
   fall in (ties: the lowest class). This affects 12 reactions of brenda_merged and
   3 of e_coli_core; the notes key `preferred_ec` of brenda_merged is not used, because a general
   SBML file does not carry it.
10. **No SBO term on reactions (R-HSA-77289).** The declared fallback is the EC class (R5). 43 of its 47
    reactions have an EC number; the four without stay outside the type boxes or in an undivided compartment.
    In e_coli_core the SBO term does not divide the cytosol (49 biochemical reactions, 1 biomass
    production), so the cytosol is also divided by EC class; its biomass reaction sits directly in the cytosol.
11. **Single-compartment model.** In `out/brenda_merged` level 1 is the EC class, not localisation, and
    the eight nodes no EC box takes are in a level-1 group named "default" of basis localisation. The
    result there is a pure type hierarchy, the same as a type-first strategy would give.
12. **Thresholds are justified, not derived from the data** (MAX_BOX 50, MIN_BOX 5, MIN_RX 3,
    HUB_BOXES 3). `RULES.md` shows how the panel totals move when each is changed.
13. **`--sbml` is not read.** Everything the rules need is in the feature tables. The file is accepted
    because the contract's command line has it.
14. **One embedded vocabulary.** The names of the seven top-level EC classes are written in
    `strategy.py` (a constant of the EC nomenclature, used only to print "EC 1 oxidoreductases" instead
    of "EC 1"). The identifier behind each label is still the EC class (`name_xref` = `ec-code:1`).
