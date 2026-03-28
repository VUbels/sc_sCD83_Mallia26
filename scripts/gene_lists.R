broad_gene_list = c("PDGFRA", "FGF7", "VIM", "EDN3", "WNT5A", "KRT1", "KRT10", "KRT6A",
                    "KRT17", "KRT75", "KRT35", "KRT85", "PECAM1", "VWF", "TAGLN", "DSP",
                    "MLANA", "PMEL", "SOX7","SOX9", "SOX18","SERPINA3", "PDZRN3", "FGL2",
                    "CXCL14", "LGR5", "COMP", "CD34", "AQP3", "KRT79", "KRT19", "WNT3",
                    "WNT10A", "WNT10B", "RGS5", "CD3D", "CD69", "TWIST2", "SOX10", "CDH9",
                    "ACTA2", "EGF", "LGR6", "FGF14", "FMN2", "CTNNA2", "GPR183", "CD53",
                    "CD2", "HLA-DQA2", "S100A2", "NES", "RGS5")

gene_list <- list(
  keratinocytes = list(
    basal = c("KRT14", "KRT5", "KRT15", "COL17A1"),
    HFSCs = c("SOX9", "LHX2", "NFATC1", "TCF3", "ITGA6", "ITGB1", "CD200", "LGR5", "LGR6", "FRZB", "FZD1", "FZD5", "FZD10", "IL31RA", "IL13RA1", "OSMR", "CD34", "CDH3", "RUNX1"),
    hair_germ = c("CD34", "CDH3", "LGR5", "CDKN2A", "RUNX1"),
    matrix = c("KRT81", "KRT83", "HOXC13", "LEF1"),
    sheath = c("KRT71", "KRT75"),
    spinous = c("KRT1", "KRT10"),
    granular = c("DSC1", "KRT2", "IVL", "TGM3"),
    RUNX3_high = c("RUNX1", "RUNX2", "RUNX3", "KRT23", "KRT18"),
    HF = c("KRT81", "KRT82", "LEF1", "KRT75", "SOX9", "HOXC13", "MGST1", "COL14A1", "CD200", "ELOVL3", "KRT79"),
    glandular = c("SCGB2A2", "SCGB1D2", "KRT7", "KRT8", "KRT19", "AQP5", "KRT79"),
    infundibulum = c("LRIG1", "KRT79", "KRT17", "PLET1", "SPRR1A", "CST6", "EFNB2", "SCA1"),
    isthmus = c("LGR6", "LRIG1", "PLET1", "KRT79", "KRT17", "PRDM1", "KRT6A", "FOXC1"),
    proliferating = c("MKI67", "CDK1", "TOP2A")
  ),
  
  endothelial = list(
    general = c("VWF", "PECAM1", "SELE", "FLT1"),
    lymphatic = c("ACKR2", "FLT4", "LYVE1", "CCL21", "TFF3", "APOD"),
    angiogenic = c("SEMA3G", "FBLN5", "NEBL", "CXCL12", "PCSK5", "SPECC1")
  ),
  
  lymphoid = list(
    Tregs = c("FOXP3", "CD4", "IL2RA", "IKZF2"),
    Th1 = c("CCR1", "CCR5", "CXCR3", "TNF", "LTA", "TBX21"),
    Th17 = c("RORA", "CCL20", "BATF", "IL1RL1", "IL6R", "IL17A", "IL17F", "IL21", "IL22"),
    Th22 = c("AHR", "CCR4", "CCR6", "CCR10", "IL6R", "IL10", "IL13", "IL22"),
    TFH = c("IL21", "CXCL13", "IL17A", "IL17F"),
    memory_T = c("CD44", "IL7R", "CCR7", "BCL6"),
    CD8NK = c("CD8A", "KLRK1", "KLRD1", "IFNG", "CCL5", "GZMA", "GZMB"),
    tissue_resident = c("CCR7", "ITGAE", "SELL", "KLRG1", "CCR4", "CCR8", "CCR10", "SELPLG")
  ),
  
  myeloid = list(
    mast = c("KIT", "ENPP3", "FCER1A", "IL1RL1", "TPSB2"),
    macrophages = c("CD163", "LGMN", "FCGR2A", "C1QB", "C5AR1", "MAFB", "FOLR2"),
    M1 = c("CCL20", "CXCL3", "IL1B", "IL6", "IL12A", "IFNG", "TNF"),
    M2a = c("CD163", "CD200", "IRF4", "TGFB1", "TGFB2", "CCL2", "STAT6"),
    TREM2 = c("TREM2", "C3", "FCGBP", "FCGR3A", "OSM", "APOE"),
    langerhans = c("ITGAX", "CD1A", "CLEC1A", "CD207", "EPCAM"),
    pDC = c("CCR7", "PTPRC", "CD209", "CLEC4C"),
    moDC = c("CD14", "CD1A", "CD1C", "ITGAX", "ITGAM", "SIRPA"),
    cDC1 = c("BTLA", "ITGAE", "CD1A", "ITGAM", "CLEC9A", "XCR1", "THBD"),
    cDC2 = c("CD14", "CD163", "CLEC10A", "NOTCH2", "ITGAM", "SIRPA", "CX3CR1", "CD1C", "CD2"),
    TCR_macs = c("CD3D", "TRAC", "TRBC1", "SPOCK2", "CD14", "CD2")
  ),
  
  fibroblasts = list(
    general = c("THY1", "COL1A1", "COL1A2", "COL3A1", "DCN", "MGP", "COL6A2", "CEBPB", "APOD", "CFD"),
    HF_associated = c("APCDD1", "VCAN", "CORIN", "PTGDS", "SOX2", "COL11A1"),
    immune_recruiting = c("CXCL1", "CXCL2", "CXCL14", "CD44"),
    papillary = c("COL6A5", "APCDD1", "HSPB3", "WIF1", "ENTPD1"),
    reticular = c("CD36"),
    dermal_papilla = c("WNT5A", "BMP4", "BMP7", "HHIP", "PTCH1", "SOX18", "RUNX1", "RUNX3", "ALX4")
  )
)