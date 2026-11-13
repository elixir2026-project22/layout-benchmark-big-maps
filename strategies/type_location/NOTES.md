# type_location - notes

All numbers below were computed from the files in `out/` (default parameters: `theta` = 1/3, `n_min` = 3,
`k_max` = 11). "Nodes" are species plus reactions. A clone is a further copy of a species listed in `clone_in`.

## Overview

| model | nodes | L1 | L2 | L3 | largest L1 share | largest L2 share | unnamed groups | L1 groups of 1-2 nodes | L2 groups of 1-2 nodes | species with clones | clones | participant edges inside the home L1 box |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| brenda_merged | 396 | 7 | 28 | 0 | 0.29 | 0.10 | 0 of 35 | 0 | 1 | 120 of 178 | 188 | 574 of 914 (63%) |
| e_coli_core | 167 | 10 | 12 | 0 | 0.22 | 0.22 | 0 of 22 | 1 | 1 | 67 of 72 | 126 | 192 of 360 (53%) |
| R-HSA-77289 | 148 | 8 | 12 | 0 | 0.25 | 0.24 | 0 of 20 | 0 | 1 | 42 of 101 | 52 | 182 of 256 (71%) |
| skm_pss | 1817 | 6 | 52 | 0 | 0.28 | 0.23 | 0 of 58 | 0 | 24 | 292 of 961 | 356 | 1853 of 2257 (82%) |
| skm_ja | 199 | 5 | 24 | 0 | 0.31 | 0.15 | 0 of 29 | 0 | 8 | 26 of 117 | 32 | 183 of 218 (84%) |
| skm_sa | 352 | 5 | 25 | 0 | 0.34 | 0.24 | 0 of 30 | 0 | 8 | 61 of 193 | 74 | 316 of 404 (78%) |

All six outputs pass `python common/contract.py validate` (`valid`). Two runs of `strategy.py` (the second
with another `PYTHONHASHSEED`) give byte-identical `groups.tsv` and `groups_meta.tsv` for all six models.
`*_grouped.xml` were written with `python common/contract.py write-sbml` and read back with libsbml: the
number of groups matches. `R-HSA-77289_grouped.xml` carries one libsbml error (a `<layout>` without
`layout:id`) that is already in the input model.

## Per model

### brenda_merged

178 species, 218 reactions. Groups: 7 at L1, 28 at L2, none at L3. Unnamed groups: 0 of 35 (0%); name sources: ec 34, sbo 1.

- L1 sizes in nodes: min 8, median 47, max 114 (29% of the model). By size class - 1: 0, 2: 0, 3-5: 0, 6-20: 2, 21-100: 3, >100: 2.
- L2 sizes in nodes: min 2, median 10.5, max 40 (10% of the model). By size class - 1: 0, 2: 1, 3-5: 4, 6-20: 15, 21-100: 8, >100: 0. 8 nodes sit directly in an L1 group without L2 (R15).
- Rules that placed the nodes: R1 4, R8 214, R11 178.
- Boundary nodes: 0, of which in an `@boundary` group: 0.
- Clones: 120 of 178 species have at least one copy; 188 copies in all, at most 5 for one species; 7 copies point to an L1 group because the type has no group for the species' compartment. Participant edges whose species is at home in the reaction's L1 box: 574 of 914; with the clones: 914 of 914.

L1 groups (nodes / reactions):

- EC 1 Oxidoreductases (`T_EC_1`, `ec-code:1`): 114 / 57
- EC 2 Transferases (`T_EC_2`, `ec-code:2`): 110 / 60
- EC 3 Hydrolases (`T_EC_3`, `ec-code:3`): 55 / 35
- EC 4 Lyases (`T_EC_4`, `ec-code:4`): 47 / 28
- EC 5 Isomerases (`T_EC_5`, `ec-code:5`): 20 / 15
- EC 6 Ligases (`T_EC_6`, `ec-code:6`): 42 / 19
- other biochemical reaction (`T_SBO_0000176_other`, `SBO:0000176`): 8 / 4

Three example groups:

- **EC 2 Transferases / EC 2.2 Transferring aldehyde or ketonic groups** (`T_EC_2__EC_2_2`, 3 reactions, 3 species; named from `ec-code:2.2`)
  - reaction `R_rx_058ab7bb35`: BR43957
  - reaction `R_rx_39f03d4644`: BR43958
  - reaction `R_rx_c3d38d6dfa`: BS320170
  - species: D-erythrose 4-phosphate (`M_mt_275a659df3`); D-xylulose 5-phosphate (`M_mt_469c9ab974`); sedoheptulose 7-phosphate (`M_mt_b8d188b145`)
- **EC 1 Oxidoreductases / EC 1.3 Acting on the CH-CH group of donors** (`T_EC_1__EC_1_3`, 4 reactions, 5 species; named from `ec-code:1.3`)
  - reaction `R_rx_08b3253feb`: BR50094
  - reaction `R_rx_0e10fa96aa`: BS370908
  - reaction `R_rx_4a2da21816`: BR49374
  - reaction `R_rx_dee9b5c8f5`: BS413615
  - species: FADH2 (`M_mt_4703982017`); crotonyl-CoA (`M_mt_4db5d73987`); FAD (`M_mt_66fd1a2692`); oxidized acceptor (`M_mt_NAME_00000036`); butyryl-CoA (`M_mt_cc8291cf5b`)
- **EC 4 Lyases / EC 4.3 Carbon-nitrogen lyases** (`T_EC_4__EC_4_3`, 7 reactions, 4 species; named from `ec-code:4.3`)
  - reaction `R_rx_56b84e5a3a`: BR3431
  - reaction `R_rx_7761ebea74`: BR481
  - reaction `R_rx_8963fc40e6`: BR3426
  - reaction `R_rx_b58b3ea39e`: BR12741
  - reaction `R_rx_c9332e8eb5`: BR487
  - reaction `R_rx_d79f5bb30e`: BR3432
  - ... and 1 more reactions
  - species: 5-aminoimidazole-4-carboxamide ribonucleotide (`M_mt_3f98b736e2`); fumarate (`M_mt_6ea7392d87`); L-threo-3-methylaspartate (`M_mt_b5f144d688`); L-arginine (`M_mt_f25ff329c6`)

### e_coli_core

72 species, 95 reactions. Groups: 10 at L1, 12 at L2, none at L3. Unnamed groups: 0 of 22 (0%); name sources: compartment 12, ec 6, sbo 4.

- L1 sizes in nodes: min 1, median 11.5, max 37 (22% of the model). By size class - 1: 1, 2: 0, 3-5: 2, 6-20: 3, 21-100: 4, >100: 0.
- L2 sizes in nodes: min 1, median 7.5, max 37 (22% of the model). By size class - 1: 1, 2: 0, 3-5: 3, 6-20: 4, 21-100: 4, >100: 0.
- Rules that placed the nodes: R1 48, R8 47, R11 72.
- Boundary nodes: 25, of which in an `@boundary` group: 25.
- Clones: 67 of 72 species have at least one copy; 126 copies in all, at most 6 for one species; 0 copies point to an L1 group because the type has no group for the species' compartment. Participant edges whose species is at home in the reaction's L1 box: 192 of 360; with the clones: 360 of 360.

L1 groups (nodes / reactions):

- EC 1 Oxidoreductases (`T_EC_1`, `ec-code:1`): 28 / 15
- EC 2 Transferases (`T_EC_2`, `ec-code:2`): 27 / 13
- EC 3 Hydrolases (`T_EC_3`, `ec-code:3`): 4 / 3
- EC 4 Lyases (`T_EC_4`, `ec-code:4`): 16 / 8
- EC 5 Isomerases (`T_EC_5`, `ec-code:5`): 7 / 5
- EC 6 Ligases (`T_EC_6`, `ec-code:6`): 5 / 3
- other biochemical reaction (`T_SBO_0000176_other`, `SBO:0000176`): 6 / 2
- translocation reaction (`T_SBO_0000185`, `SBO:0000185`): 36 / 25
- exchange reaction (`T_SBO_0000627`, `SBO:0000627`): 37 / 20
- biomass production (`T_SBO_0000629`, `SBO:0000629`): 1 / 1

Three example groups:

- **EC 5 Isomerases / cytosol** (`T_EC_5__c`, 5 reactions, 2 species; named from `c`, compartment_id `c`)
  - reaction `R_PGI`: Glucose-6-phosphate isomerase
  - reaction `R_PGM`: Phosphoglycerate mutase
  - reaction `R_RPE`: Ribulose 5-phosphate 3-epimerase
  - reaction `R_RPI`: Ribose-5-phosphate isomerase
  - reaction `R_TPI`: Triose-phosphate isomerase
  - species: 3-Phospho-D-glycerate (`M_3pg_c`); D-Ribulose 5-phosphate (`M_ru5p__D_c`)
- **EC 6 Ligases / cytosol** (`T_EC_6__c`, 3 reactions, 2 species; named from `c`, compartment_id `c`)
  - reaction `R_GLNS`: Glutamine synthetase
  - reaction `R_GLUN`: Glutaminase
  - reaction `R_SUCOAS`: Succinyl-CoA synthetase (ADP-forming)
  - species: L-Glutamine (`M_gln__L_c`); Ammonium (`M_nh4_c`)
- **translocation reaction / cytosol boundary** (`T_SBO_0000185__c__boundary`, 25 reactions, 0 species; named from `c`, compartment_id `c@boundary`)
  - reaction `R_ACALDt`: Acetaldehyde reversible transport
  - reaction `R_ACt2r`: Acetate reversible transport via proton symport
  - reaction `R_AKGt2r`: 2 oxoglutarate reversible transport via symport
  - reaction `R_ATPS4r`: ATP synthase (four protons for one ATP)
  - reaction `R_CO2t`: CO2 transporter via diffusion
  - reaction `R_CYTBD`: Cytochrome oxidase bd (ubiquinol-8: 2 protons)
  - ... and 19 more reactions

### R-HSA-77289

101 species, 47 reactions. Groups: 8 at L1, 12 at L2, none at L3. Unnamed groups: 0 of 20 (0%); name sources: compartment 12, ec 7, sbo 1.

- L1 sizes in nodes: min 3, median 17, max 37 (25% of the model). By size class - 1: 0, 2: 0, 3-5: 1, 6-20: 4, 21-100: 3, >100: 0.
- L2 sizes in nodes: min 1, median 11, max 36 (24% of the model). By size class - 1: 1, 2: 0, 3-5: 3, 6-20: 7, 21-100: 1, >100: 0.
- Rules that placed the nodes: R6 2, R8 41, R10 4, R11 101.
- Boundary nodes: 4, of which in an `@boundary` group: 4.
- Clones: 42 of 101 species have at least one copy; 52 copies in all, at most 4 for one species; 0 copies point to an L1 group because the type has no group for the species' compartment. Participant edges whose species is at home in the reaction's L1 box: 182 of 256; with the clones: 256 of 256.

L1 groups (nodes / reactions):

- EC 1.1 Acting on the CH-OH group of donors (`T_EC_1_1`, `ec-code:1.1`): 17 / 7
- EC 1.3 Acting on the CH-CH group of donors (`T_EC_1_3`, `ec-code:1.3`): 37 / 11
- EC 1.6 Acting on NADH or NADPH (`T_EC_1_6`, `ec-code:1.6`): 3 / 1
- EC 2.3 Acyltransferases (`T_EC_2_3`, `ec-code:2.3`): 23 / 8
- EC 4.2 Carbon-oxygen lyases (`T_EC_4_2`, `ec-code:4.2`): 16 / 7
- EC 5 Isomerases (`T_EC_5`, `ec-code:5`): 10 / 4
- EC 6 Ligases (`T_EC_6`, `ec-code:6`): 17 / 3
- other biochemical reaction (`T_SBO_0000176_other`, `SBO:0000176`): 25 / 6

Three example groups:

- **EC 5 Isomerases / mitochondrial matrix** (`T_EC_5__compartment_5460`, 4 reactions, 6 species; named from `compartment_5460`, compartment_id `compartment_5460`)
  - reaction `reaction_109338`: Isomerization of cis,cis-3,6-Dodecadienoyl-CoA to form trans,cis-Lauro-2,6-dienoyl-CoA
  - reaction `reaction_109998`: Isomerization of 3-trans-decenoyl-CoA to form trans-dec-2-enoyl-CoA
  - reaction `reaction_71010`: MUT isomerises L-MM-CoA to SUCC-CoA
  - reaction `reaction_71020`: D-methylmalonyl-CoA <=> L-methylmalonyl-CoA
  - species: DCI dimer [mitochondrial matrix] (`species_110060`); 2xMMAA:2xMUT:AdoCbl [mitochondrial matrix] (`species_3159272`); SUCC-CoA [mitochondrial matrix] (`species_70958`); L-MM-CoA [mitochondrial matrix] (`species_71009`); MCEE [mitochondrial matrix] (`species_71012`); trans,cis-Lauro-2,6-dienoyl-CoA  [mitochondri (`species_77360`)
- **EC 2.3 Acyltransferases / mitochondrial inner membrane** (`T_EC_2_3__compartment_5460__boundary`, 1 reactions, 2 species; named from `compartment_70616`, compartment_id `compartment_5460@boundary`)
  - reaction `reaction_8933547`: MCAT transfers Mal from Mal-CoA to NDUFAB1
  - species: NDUFAB1 [mitochondrial inner membrane] (`species_163925`); Mal-PPANT-S112-NDUFAB1 [mitochondrial inner m (`species_8934586`)
- **other biochemical reaction / cytosol** (`T_SBO_0000176_other__compartment_70101`, 2 reactions, 10 species; named from `compartment_70101`, compartment_id `compartment_70101`)
  - reaction `reaction_5690043`: Cytosolic ACOTs hydrolyse MCFA-CoA, LCFA-CoA
  - reaction `reaction_8848250`: ACBD6 binds LCFA-CoA
  - species: LCFA-CoA [cytosol] (`species_164991`); H2O [cytosol] (`species_29356`); PCTP [cytosol] (`species_5215990`); MCFA-CoA, LCFA-CoA [cytosol] (`species_5690044`); Cytosolic ACOTs, THEM4 [cytosol] (`species_5690052`); MCFA, LCFA [cytosol] (`species_5690064`); ... and 4 more

### skm_pss

961 species, 856 reactions. Groups: 6 at L1, 52 at L2, none at L3. Unnamed groups: 0 of 58 (0%); name sources: compartment 52, sbo 6.

- L1 sizes in nodes: min 41, median 371, max 500 (28% of the model). By size class - 1: 0, 2: 0, 3-5: 0, 6-20: 0, 21-100: 1, >100: 5.
- L2 sizes in nodes: min 1, median 4, max 418 (23% of the model). By size class - 1: 12, 2: 12, 3-5: 4, 6-20: 11, 21-100: 9, >100: 4.
- Rules that placed the nodes: R1 855, R10 1, R11 961.
- Boundary nodes: 110, of which in an `@boundary` group: 110.
- Clones: 292 of 961 species have at least one copy; 356 copies in all, at most 4 for one species; 3 copies point to an L1 group because the type has no group for the species' compartment. Participant edges whose species is at home in the reaction's L1 box: 1853 of 2257; with the clones: 2257 of 2257.

L1 groups (nodes / reactions):

- other biochemical reaction (`T_SBO_0000176_other`, `SBO:0000176`): 387 / 157
- non-covalent binding (`T_SBO_0000177`, `SBO:0000177`): 500 / 186
- degradation (`T_SBO_0000179`, `SBO:0000179`): 41 / 26
- conversion (`T_SBO_0000182`, `SBO:0000182`): 404 / 198
- genetic production (`T_SBO_0000589`, `SBO:0000589`): 355 / 206
- transport reaction (`T_SBO_0000655`, `SBO:0000655`): 130 / 83

Three example groups:

- **degradation / nucleus** (`T_SBO_0000179__c_nuc`, 8 reactions, 3 species; named from `c_nuc`, compartment_id `c_nuc`)
  - reaction `rx00017`: EIN3(like) -> (none) [modifiers: EBF|SCF]
  - reaction `rx00044`: JAZ -> (none) [modifiers: COI1|JA-Ile|SCF]
  - reaction `rx00112`: ARF -> (none) [modifiers: TAS3]
  - reaction `rx00123`: WRKY53 -> (none) [modifiers: UPL5]
  - reaction `rx00150`: CAMTA3 -> (none) [modifiers: SR1IP1]
  - reaction `rx00422`: NPR1 -> (none) [modifiers: NPR3, SA]
  - ... and 2 more reactions
  - species: NPR3 (`s_NPR3_nuc_pa`); NPR4 (`s_NPR4_nuc_pa`); SR1IP1 (`s_SR1IP1_nuc_pa`)
- **conversion / apoplast** (`T_SBO_0000182__c_ap`, 5 reactions, 2 species; named from `c_ap`, compartment_id `c_ap`)
  - reaction `rx00174`: RBOH -> RBOH [modifiers: CDPK]
  - reaction `rx00239`: RBOH -> RBOH [modifiers: Ca2+]
  - reaction `rx00427`: RBOH -> RBOH [modifiers: BIK1]
  - reaction `rx00449`: RBOH -> RBOH [modifiers: LECRK19]
  - reaction `rx00451`: RBOH -> RBOH [modifiers: SIK1]
  - species: RBOH (`s_RBOH_ap_p`); RBOH (`s_RBOH_ap_pa`)
- **transport reaction / chloroplast boundary** (`T_SBO_0000655__c_chl__boundary`, 12 reactions, 0 species; named from `c_chl`, compartment_id `c_chl@boundary`)
  - reaction `rx00396`: OPDA -> OPDA [modifiers: JASSY]
  - reaction `rx00456`: EIJP1 -> EIJP1
  - reaction `rx00517`: IsoChor -> IsoChor [modifiers: EDS5]
  - reaction `transport_CA_m_chl_to_cyt`: CA -> CA
  - reaction `transport_EDS5_pa_cyt_to_chl`: EDS5 -> EDS5
  - reaction `transport_IAA_m_cyt_to_chl`: IAA -> IAA
  - ... and 6 more reactions

### skm_ja

117 species, 82 reactions. Groups: 5 at L1, 24 at L2, none at L3. Unnamed groups: 0 of 29 (0%); name sources: compartment 24, sbo 5.

- L1 sizes in nodes: min 19, median 38, max 61 (31% of the model). By size class - 1: 0, 2: 0, 3-5: 0, 6-20: 1, 21-100: 4, >100: 0.
- L2 sizes in nodes: min 1, median 4, max 30 (15% of the model). By size class - 1: 5, 2: 3, 3-5: 7, 6-20: 6, 21-100: 3, >100: 0.
- Rules that placed the nodes: R1 81, R10 1, R11 117.
- Boundary nodes: 16, of which in an `@boundary` group: 16.
- Clones: 26 of 117 species have at least one copy; 32 copies in all, at most 4 for one species; 1 copies point to an L1 group because the type has no group for the species' compartment. Participant edges whose species is at home in the reaction's L1 box: 183 of 218; with the clones: 218 of 218.

L1 groups (nodes / reactions):

- other biochemical reaction (`T_SBO_0000176_other`, `SBO:0000176`): 61 / 21
- non-covalent binding (`T_SBO_0000177`, `SBO:0000177`): 36 / 12
- conversion (`T_SBO_0000182`, `SBO:0000182`): 45 / 20
- genetic production (`T_SBO_0000589`, `SBO:0000589`): 38 / 19
- transport reaction (`T_SBO_0000655`, `SBO:0000655`): 19 / 10

Three example groups:

- **transport reaction / nucleus boundary** (`T_SBO_0000655__c_nuc__boundary`, 4 reactions, 0 species; named from `c_nuc`, compartment_id `c_nuc@boundary`)
  - reaction `transport_COI1_pa_cyt_to_nuc`: COI1 -> COI1
  - reaction `transport_JAM_pa_cyt_to_nuc`: JAM -> JAM
  - reaction `transport_JAZ_pa_nuc_to_cyt`: JAZ -> JAZ
  - reaction `transport_MYC2_pa_cyt_to_nuc`: MYC2 -> MYC2
- **non-covalent binding / nucleus** (`T_SBO_0000177__c_nuc`, 3 reactions, 9 species; named from `c_nuc`, compartment_id `c_nuc`)
  - reaction `rx00018`: EIN3(like) + JAZ -> EIN3(like)|JAZ
  - reaction `rx00159`: DELLA + JAZ -> DELLA|JAZ
  - reaction `rx00162`: MYC2 + DELLA -> DELLA|MYC
  - species: COI1 (`s_COI1_nuc_pa`); DELLA|JAZ (`s_DELLAJAZ_nuc_ca`); DELLA|MYC (`s_DELLAMYC_nuc_ca`); DELLA (`s_DELLA_nuc_pa`); EIN3(like)|JAZ (`s_EIN3likeJAZ_nuc_ca`); EIN3(like) (`s_EIN3like_nuc_pa`); ... and 3 more
- **conversion / peroxisome** (`T_SBO_0000182__c_per`, 1 reactions, 3 species; named from `c_per`, compartment_id `c_per`)
  - reaction `rx00312`: ACX -> ACX [modifiers: CAT]
  - species: ACX (`s_ACX_per_p`); ACX (`s_ACX_per_pa`); CAT (`s_CAT_per_pa`)

### skm_sa

193 species, 159 reactions. Groups: 5 at L1, 25 at L2, none at L3. Unnamed groups: 0 of 30 (0%); name sources: compartment 25, sbo 5.

- L1 sizes in nodes: min 28, median 69, max 120 (34% of the model). By size class - 1: 0, 2: 0, 3-5: 0, 6-20: 0, 21-100: 3, >100: 2.
- L2 sizes in nodes: min 1, median 6, max 83 (24% of the model). By size class - 1: 3, 2: 5, 3-5: 4, 6-20: 7, 21-100: 6, >100: 0.
- Rules that placed the nodes: R1 156, R10 3, R11 193.
- Boundary nodes: 27, of which in an `@boundary` group: 27.
- Clones: 61 of 193 species have at least one copy; 74 copies in all, at most 3 for one species; 0 copies point to an L1 group because the type has no group for the species' compartment. Participant edges whose species is at home in the reaction's L1 box: 316 of 404; with the clones: 404 of 404.

L1 groups (nodes / reactions):

- other biochemical reaction (`T_SBO_0000176_other`, `SBO:0000176`): 105 / 38
- non-covalent binding (`T_SBO_0000177`, `SBO:0000177`): 69 / 24
- conversion (`T_SBO_0000182`, `SBO:0000182`): 28 / 13
- genetic production (`T_SBO_0000589`, `SBO:0000589`): 120 / 65
- transport reaction (`T_SBO_0000655`, `SBO:0000655`): 30 / 19

Three example groups:

- **other biochemical reaction / chloroplast** (`T_SBO_0000176_other__c_chl`, 4 reactions, 8 species; named from `c_chl`, compartment_id `c_chl`)
  - reaction `rx00069`: Chor -> IsoChor [modifiers: ICS]
  - reaction `rx00073`: Chor -> Prep [modifiers: CM]
  - reaction `rx00076`: Phe -> CA [modifiers: PAL]
  - reaction `rx00541`: L-arogenate -> Phe [modifiers: ADT]
  - species: ADT (`s_ADT_chl_pa`); CA (`s_CA_chl_m`); CM (`s_CM_chl_pa`); Chor (`s_Chor_chl_m`); IsoChor (`s_IsoChor_chl_m`); PAL (`s_PAL_chl_pa`); ... and 2 more
- **transport reaction / chloroplast boundary** (`T_SBO_0000655__c_chl__boundary`, 6 reactions, 0 species; named from `c_chl`, compartment_id `c_chl@boundary`)
  - reaction `rx00517`: IsoChor -> IsoChor [modifiers: EDS5]
  - reaction `transport_CA_m_chl_to_cyt`: CA -> CA
  - reaction `transport_EDS5_pa_cyt_to_chl`: EDS5 -> EDS5
  - reaction `transport_Larogenate_m_cyt_to_chl`: L-arogenate -> L-arogenate
  - reaction `transport_Phe_m_chl_to_cyt`: Phe -> Phe
  - reaction `transport_Prep_m_chl_to_cyt`: Prep -> Prep
- **conversion / nucleus** (`T_SBO_0000182__c_nuc`, 5 reactions, 6 species; named from `c_nuc`, compartment_id `c_nuc`)
  - reaction `activation_MYC2_pa_nuc`: MYC2 -> MYC2
  - reaction `rx00432`: MYC2 -> MYC2 [modifiers: NPR1]
  - reaction `rx00478`: bZIP21 -> bZIP21 [modifiers: NPR1]
  - reaction `rx00511`: NPR1 -> NPR1 [modifiers: TCP8,14,15]
  - reaction `rx00544`: NPR1|TGA -> NPR1|TGA [modifiers: EDS1|PAD4|SAG]
  - species: MYC2 (`s_MYC2_nuc_p`); MYC2 (`s_MYC2_nuc_pa`); NPR1 (`s_NPR1_nuc_p`); TCP8,14,15 (`s_TCP81415_nuc_pa`); bZIP21 (`s_bZIP21_nuc_p`); bZIP21 (`s_bZIP21_nuc_pa`)

## Sensitivity to the parameters

Same program on all six models; L1 groups / L2 groups / largest L1 share of nodes.

| model | default | theta = 1/2 | n_min = 2 | n_min = 5 |
|---|---|---|---|---|
| brenda_merged | 7 / 28 / 0.29 | 7 / 28 / 0.29 | 7 / 29 / 0.29 | 7 / 23 / 0.29 |
| e_coli_core | 10 / 12 / 0.22 | 10 / 12 / 0.22 | 10 / 12 / 0.22 | 8 / 10 / 0.22 |
| R-HSA-77289 | 8 / 12 / 0.25 | 6 / 10 / 0.39 | 10 / 14 / 0.25 | 6 / 10 / 0.35 |
| skm_pss | 6 / 52 / 0.28 | 6 / 52 / 0.28 | 6 / 52 / 0.28 | 6 / 52 / 0.28 |
| skm_ja | 5 / 24 / 0.31 | 5 / 24 / 0.31 | 5 / 24 / 0.31 | 5 / 24 / 0.31 |
| skm_sa | 5 / 25 / 0.34 | 3 / 20 / 0.58 | 6 / 26 / 0.34 | 5 / 25 / 0.34 |

`k_max` = 11 never binds: the outputs with `k_max` = 99 are byte-identical.

## Known failure cases

1. **A reaction type is not a pathway.** Consecutive steps of a metabolic route usually belong to different
   enzyme classes, so the route is cut at almost every step. In e_coli_core glycolysis is spread over
   isomerases, transferases, lyases and oxidoreductases; only 192 of 360 participant edges
   (53%) stay inside the home L1 box, in brenda_merged 574 of 914 (63%). Agreement
   with curated pathway labels should be expected to be low for the metabolic models; this is a property of
   the hierarchy, not of a parameter.
2. **Many clones in metabolic models.** 67 of 72 species of e_coli_core and 120 of 178 of brenda_merged
   have at least one copy. Currency metabolites are copied into nearly every type box (copies in
   e_coli_core: H2O (6), H+ (6), Phosphate (6), ADP (5), ATP (5)). The signalling models need fewer (26 of 117 in skm_ja,
   61 of 193 in skm_sa, 292 of 961 in skm_pss).
3. **Tiny localisation groups.** Type x space leaves sparse cells: 24 of 52 L2 groups of skm_pss hold one
   or two nodes (9 species only, 11 reactions only, 4 mixed; 8 of them are `@boundary` groups), 8 of 24
   in skm_ja, 8 of 25 in skm_sa. They cannot be merged without breaking the localisation rule: a vacuole
   with one species, or a boundary crossed by one reaction of a type, is a group of one.
4. **Tiny L1 groups.** e_coli_core: "biomass production" is one reaction and no species (all its participants
   are at home elsewhere). R-HSA-77289: "EC 1.6 Acting on NADH or NADPH" is one reaction with two species,
   what was left of the oxidoreductases after EC 1.1 and EC 1.3 were split off. e_coli_core: "other biochemical
   reaction" holds the two reactions without EC number (6 nodes), "EC 3 Hydrolases" 4 nodes, "EC 6 Ligases" 5.
5. **Large groups that the strategy cannot divide.** skm_pss "non-covalent binding / cytoplasm" holds 418
   nodes (23% of the model): the type is a leaf of SBO and the location is given, so nothing is left to
   split it by. An indivisible type can also exceed `theta`: "genetic production" is 65 of 159 reactions
   (41%) of skm_sa.
6. **"other biochemical reaction" is named but weak.** It is what SBO and EC do not specify further: 157
   reactions and 387 nodes (21%) of skm_pss, mostly catalysed metabolic steps annotated with the
   plain term. A reader can say what the box is ("biochemical reactions without a finer type") but not what
   it does.
7. **"conversion" is not split** into activation and inactivation, and "transport reaction" not into active
   and passive, because neither exceeds `theta`. A reader who wants the signalling sign at the top level
   does not get it.
8. **Reactions with several EC numbers.** 113 of 214 EC-annotated reactions of brenda_merged carry two or more EC
   numbers, 12 of them from more than one class (15 and 3 of 51 in e_coli_core). The majority vote of R8
   can misplace them: e_coli_core "Glutaminase" (`R_GLUN`, a hydrolase) carries nine EC numbers, four of them
   in class 6, and lands in "EC 6 Ligases".
9. **The fallback for reactions without SBO term is an inference.** R-HSA-77289 has no reaction SBO term:
   41 reactions are typed by their EC class (R8), two by R6 ("biochemical reaction", lumped Reactome steps
   without EC), and four are pooled by R10 (two hydrolases, two binding steps recognised by R4) into "other
   biochemical reaction", which therefore mixes thioesterase steps, two lumped beta-oxidation steps and two
   binding steps. R6 would also call any single-space reaction with reactants and products a biochemical
   reaction. The split by EC subclass works on this pathway because its steps differ in enzyme class; a
   pathway whose steps share one class would stay one box.
10. **Fallback rules not exercised by the panel.** R2, R3, R5, R7 and R12 placed no node of the six models.
    They were run only on a nine-node synthetic model (valid output), so their behaviour on real data is untested.
11. **One compartment only (brenda_merged).** L2 is not localisation there but the EC subclass (R15), which
    departs from "L2 = localisation"; the alternative was one L2 group named "default" inside every type.
    The 8 nodes of "other biochemical reaction" (no EC number) have no L2. Two L2 groups there hold a single
    reaction with its species ("EC 4.6", "EC 6.1").
12. **Clones without a localisation group.** 7 copies in brenda_merged, 3 in skm_pss and 1 in skm_ja point to
    an L1 group, because that type has no group for the species' own compartment (in brenda_merged: the L1
    group without L2).
13. **Small models can stay one box.** If no subtype reaches `n_min` reactions the type is not divided; a
    model of four reactions of four types is one group named "process".
14. **Repeated labels.** The label of a localisation group is the compartment name, so "cytoplasm" appears
    once per type; the groups are told apart by id and, in the SBML, by the path name "type / compartment".
15. **Names depend on two embedded tables.** EC class and subclass names and six SBO names for the fallback
    are constants in `strategy.py`, not read from the model. The `--sbml` file is not read at all.
