## Manually input Qiagen list miR21 targets mouse
Targets <- c("Basp1", "Bmpr2", "Btg2", "Cdc25a", "Cdk2ap1", "Daxx", "Derl1", "E2f1", "E2f2", "Eif2s1",
             "Eif4a2", "Fasl", "Hnrnpk", "Icam1", "Il1b", "Jag1", "Jmy", "Lrrfip1", "Marcks", "Mef2c",
             "Msh2", "Msh6", "Mtap", "Myc", "Ncapg", "Ncoa3", "Pcbp1", "Pdcd4", "Pdha2", "Peli1", "Plat",
             "Plod3", "Ppif", "Ppp1r3b", "Pten", "Rasgrp1", "Reck", "Rhob", "Rps7", "Rtn4", "Serpinb5",
             "Sox5", "Spats2l", "Spry1", "Spry2", "Tgfbi", "Tgfbr2", "Tgfbr3", "Tiam1", "Timp3", "Tm9sf3",
             "Tnrc6b", "Topors", "Tpm1", "Trp53bp2", "Trp63", "Wfs1", "Wibg", "Acat1",  "Acvr2a", "Adnp",
             "Ankrd46", "Apaf1", "Arhgap24", "Armcx1", "Asf1a", "Bcl11b", "Bcl2", "Cep68", "Crebl2",
             "Dag1", "Egr3", "Elf2", "Klf5", "Klhl42", "Krit1", "Lancl1", "Nfib", "Ntf3", "Pitx2", "Ppp1r3a",
             "Rasa1", "Tnks", "Ube2d3")

# Mir21 miRTarBase Strong Evidence
miRTarBase_miR21_strong = c("Pdcd4", "Pten", "Reck", "Tnfaip8l2", "Fasl", "Peli1", "Spry2", "Spry1", 
                            "Tgfbi", "Tgfbr3", "Smad7", "Yod1", "Sod2", "Sox5", "Rc3h1", "Spry4", "Spry3",
                            "Gt(ROSA)26Sor", "Yy1", "Eif4e3", "Pdcd10", "Mmp9", "Pias3", "Timp3", "Map2k3", 
                            "Sorbs2", "Sox2")

miRTarBase_miR21_strong_order_CP4_TP5 = c("Yod1", "Smad7", "Pdcd4", "Pten", "Tnfaip8l2", "Yy1", "Pias3", "Rc3h1", 
                                          "Tgfbr3", "Peli1", "Pdcd10", "Spry3", "Map2k3", "Sox5", "Fasl", "Spry1",
                                          "Mmp9", "Eif4e3", "Gt(ROSA)26Sor", "Reck", "Sorbs2", "Timp3", "Spry2", 
                                          "Sod2", "Tgfbi", "Spry4", "Sox2")


# Mir21 miRTarBase Strong Evidence (Excel file available from miRTarBase supplemented with website (miRTarBase) search results). Gene removed without evidence.
miRTarBase_miR21_strong_checked = c("Fasl", "Peli1","Pdcd4","Spry2","Pten","Reck","Spry1","Tgfbi","Pias3","Tgfbr3","Tnfaip8l2","Smad7",
                            "Gt(ROSA)26Sor","Yy1","Eif4e3","Pdcd10","Timp3","Yod1","Mmp9","Map2k3","Sod2","Sorbs2","Sox5","Rc3h1", "Sox2")

# Mir21 miRTarBase Strong Evidence + Own Validated Targets
miRTarBase_miR21_strong_extended = c("Pdcd4", "Pten", "Reck", "Tnfaip8l2", "Fasl", "Peli1", "Spry2", "Spry1", 
                                      "Tgfbi", "Tgfbr3", "Smad7", "Yod1", "Sod2", "Sox5", "Rc3h1", "Spry4", "Spry3",
                                      "Gt(ROSA)26Sor", "Yy1", "Eif4e3", "Pdcd10", "Mmp9", "Pias3", "Timp3", "Map2k3", 
                                      "Sorbs2", "Sox2", "Sowahc", "Galnt12", "Btg2", "Mef2c")

# Mir21 miRTarBase Strong Evidence + Own Validated Targets + Targets found in literature - Target w/o evidence
miRTarBase_miR21_strong_validated = c("Pdcd4", "Pten", "Reck", "Tnfaip8l2", "Fasl", "Peli1", "Spry2", "Spry1", 
                                      "Tgfbi", "Tgfbr3", "Smad7", "Yod1", "Sod2", "Sox5",
                                     "Yy1", "Eif4e3", "Pdcd10", "Mmp9", "Timp3", "Map2k3", 
                                      "Sox2", "Sowahc", "Galnt12", "Btg2", "Mef2c", "Bmpr2")

# Mir21 miRTarBase + miRWalk2.0 Strong Evidence + Own Validated Targets + Targets found in literature - Target w/o evidence
Targets_miR21_strong_validated = c("Yod1","Smad7","Pdcd4","Pten","Sowahc","Tnfaip8l2","Galnt12","Yy1","Tgfbr3",
                                      "Peli1","Pdcd10","Map2k3","Sox5","Fasl","Spry1","Mmp9","Eif4e3","Reck","Timp3",
                                      "Spry2","Sod2","Tgfbi","Sox2","Btg2","Bmpr2","Il12a","Rhob","Bcl7a","Reck",
                                      "Ski","Apaf1","Nfib","Marcks","Thrb","Wwp1","Sox7","Gramd3","Cntfr","Cdc25a",
                                      "Fnip1","Klf3","Cd44","Atpaf1","Dnajc16","Nfat5","Pcsk6","Lemd3","Lrrc57","Ube2d3",
                                      "Rbpj","Rmnd5a","Kbtbd2","Adgrg2","Dusp8","Prpf4b","Krit1","Cxcl10","Zadh2")

# Mir21 miRTarBase + miRWalk2.0 Strong Evidence + Own Validated Targets + Targets found in literature - Target w/o evidence
Targets_miR21_strong_validated_181011 = c("Yod1","Smad7","Pdcd4","Pten","Tnfaip8l2","Yy1","Tgfbr3","Peli1","Pdcd10","Map2k3",
                                        "Sox5","Fasl","Spry1","Mmp9","Eif4e3","Reck","Timp3","Spry2","Sod2","Tgfbi","Sox2",
                                        "Btg2","Bmpr2","Il12a","Rhob","Bcl7a","Reck","Ski","Apaf1","Nfib","Marcks","Thrb",
                                        "Wwp1","Sox7","Gramd3","Cntfr","Cdc25a","Fnip1","Klf3","Cd44","Atpaf1","Dnajc16","Nfat5",
                                        "Pcsk6","Lemd3","Lrrc57","Ube2d3","Rbpj","Rmnd5a","Kbtbd2","Adgrg2","Dusp8","Prpf4b",
                                        "Krit1","Cxcl10","Zadh2","Kbtbd7","Stat3","Cdk6","Gsk3b","Cdk2ap1","Tpm1","Bcl2")



# Mir21 miRTarBase Strong Evidence Human converted to Mouse (Use Jaxgenomics tool)

miRTarBase_miR21_strong_hsa_to_mmu = c("Akt2","Ankrd46","Anp32a","Apaf1","Basp1","Bcl2","Bcl6","Bmpr2","Btg2","Ccl20",
                                        "Ccr1","Cdc25a","Cdk2ap1","Clu","Col4a1","Daxx","Derl1","Dock4","Dock5","Dock7",
                                        "Dusp10","E2f1","Egfr","Eif4a2","Erbb2","Fasl","Fmod","Gas5","Gdf5","Hnrnpk",
                                        "Hpgd","Icam1","Igf1r","Il1b","Il1rn","Irak1","Iscu","Jag1","Jmy","Lrrfip1",
                                        "Map2k3","Marcks","Mef2c","Msh2","Msh6","Mtap","Myd88","Ncapg","Ncoa3","Nfia",
                                        "Nfib","Ntf3","Pcbp1","Pdcd4","Pias3","Plat","Plod2","Plod3","Ppara","Ppif",
                                        "Pten","Ptx3","Pitx3","Rasa1","Syngap1","Rasgrp1","Reck","Rest","Ea2","Rho",
                                        "Rhod","Rhob","Rps7","Rtn4","Satb1","Serpinb5","Serpini1","Setd2","Sirt2","Smad7",
                                        "Smarca4","Smn1","Sod3","Sox5","Sp1","Dand5","Spry2","Stat3","Tcf21","Tgfbi",
                                        "Tgfbr2","Tgfbr3","Tgif1","Tiam1","Timp3","Tm9sf3","Tnfaip3","Tnfrsf10b","Topors",
                                        "Trp53bp2","Tcp1","Trp63","Tpm1","Tcpm1","Vegfa","Vhl","Wwp1","Yod1")

# Mir21 miRTarBase Strong Evidence Mouse and Human converted to Mouse (Use Jaxgenomics tool)

miRTarBase_miR21_strong_mmu_and_hsa = c("Akt2","Ankrd46","Anp32a","Apaf1","Basp1","Bcl2","Bcl6","Bmpr2","Btg2","Ccl20",
                                        "Ccr1","Cdc25a","Cdk2ap1","Clu","Col4a1","Daxx","Derl1","Dock4","Dock5","Dock7",
                                        "Dusp10","E2f1","Egfr","Eif4a2","Erbb2","Fmod","Gas5","Gdf5","Hnrnpk",
                                        "Hpgd","Icam1","Igf1r","Il1b","Il1rn","Irak1","Iscu","Jag1","Jmy","Lrrfip1",
                                        "Marcks","Mef2c","Msh2","Msh6","Mtap","Myd88","Ncapg","Ncoa3","Nfia",
                                        "Nfib","Ntf3","Pcbp1","Plat","Plod2","Plod3","Ppara","Ppif",
                                        "Ptx3","Pitx3","Rasa1","Syngap1","Rasgrp1","Rest","Ea2","Rho",
                                        "Rhod","Rhob","Rps7","Rtn4","Satb1","Serpinb5","Serpini1","Setd2","Sirt2",
                                        "Smarca4","Smn1","Sod3","Sp1","Dand5","Stat3","Tcf21",
                                        "Tgfbr2","Tgif1","Tiam1","Tm9sf3","Tnfaip3","Tnfrsf10b","Topors",
                                        "Trp53bp2","Tcp1","Trp63","Tpm1","Tcpm1","Vegfa","Vhl","Wwp1",
                                        "Pdcd4", "Pten", "Reck", "Tnfaip8l2", "Fasl", "Peli1", "Spry2", "Spry1", 
                                      "Tgfbi", "Tgfbr3", "Smad7", "Yod1", "Sod2", "Sox5", "Rc3h1", "Spry4", "Spry3",
                                      "Gt(ROSA)26Sor", "Yy1", "Eif4e3", "Pdcd10", "Mmp9", "Pias3", "Timp3", "Map2k3", 
                                      "Sorbs2", "Sox2", "Sowahc", "Galnt12")



