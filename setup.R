# SET UP

library(hemibrainr)
library(natverse)
library(fafbseg)
library(reticulate)
library(googlesheets4)
library(dplyr)
library(neuprintr)
# library(catmaid)
# library(devtools)
# library(plotly)
# library(stringr)
library(ggplot2)
# library(nat.jrcbrains)
# library(rgl)
# library(forcats)
# library(stats)

# RGL Usermat
rgl_usermat = structure(c(0.98662239,  -0.1360703,  0.08978364, 0, 
                          0.06441419,  -0.1805410, -0.98145604, 0, 
                          0.14975676,   0.9741101, -0.16936085, 0, 
                          0,            0,          0,          1), .Dim = c(4L, 4L))
rgl_zoom = 1
rgl_rect = c(0L, 45L, 1171L, 823L)

# flywire setup
choose_segmentation("flywire31")
catmaid_login() #this only works if you have your credentials saved in your environment

# load flywiretypes
flywire_types <- readRDS(file ="/Users/tomke/Documents/dev/LH_MB_Downstream/flywire_types_05.10.Rds")

# OR update them: get all flywire neurons currently matched/typed
# MBONs
gs4_auth("ts587@cam.ac.uk")
MBON = read_sheet("1BfYo4AspojLNMAvoE5VjGlpO3Vxu8R95B21wHfKf0Mg")
MBONpts= subset(MBON$flywire.xyz,(!is.na(MBON$flywire.xyz)))
MBONid = flywire_xyz2id(MBONpts, rawcoords=TRUE)
MBONtype = subset(MBON$Type,(MBON$flywire.xyz!= 'NA'))
MBONname= data.frame(MBONtype, MBONid, MBONpts)
names(MBONname)[1] = "type"
names(MBONname)[2] = "post_id"
names(MBONname)[3] = "xyz"
# LHCENTs
LHCENT = read_sheet("1prxk5N4gekXEtGy2UdYKLCeF2XNU2I9IZr-2OJQTVjQ", sheet = 'master')
pts_R= subset(LHCENT$`Flywire R XYZ`,(!is.na(LHCENT$`Flywire R XYZ`)))
pts_L= subset(LHCENT$`Flywire L XYZ`,(!is.na(LHCENT$`Flywire L XYZ`)))
LHCENTpts = c(as.character(pts_R), as.character(pts_L))
LHCENTpts = LHCENTpts[LHCENTpts != "*"]
LHCENTid = flywire_xyz2id(LHCENTpts, rawcoords=TRUE)
LHCENTtype_R = subset(LHCENT$`FIB R name`,(!is.na(LHCENT$`Flywire R XYZ`)))
LHCENTtype_R = subset(LHCENT$`FIB R name`,(LHCENT$`Flywire R XYZ`!= "*"))
LHCENTtype_L = subset(LHCENT$`FIB R name`,(!is.na(LHCENT$`Flywire L XYZ`)))
LHCENTtype = c(as.character(LHCENTtype_R), as.character(LHCENTtype_L))
# LHCENTtype= unlist(LHCENTtype, use.names=FALSE)
LHCENTname= data.frame(LHCENTtype, LHCENTid, LHCENTpts)
names(LHCENTname)[1] = "type"
names(LHCENTname)[2] = "post_id"
names(LHCENTname)[3] = "xyz"
# DNs
gs4_auth("tomke.stuerner@googlemail.com")
putative_DN = read_sheet("1vwlsODui1_kHCvLSyEWkbOoTwBZa52nJQNV_7BVSXc4", sheet = "FlyWire_DNs_April")
DNpts= subset(putative_DN$flywire.xyz,(!is.na(putative_DN$flywire.xyz)))
# p_DNid = flywire_xyz2id(DNpts, rawcoords=TRUE)
p_DNid = putative_DN$flywire.id
DNpid = data.frame(DNpts,p_DNid)
DNtype = replicate(length(p_DNid), "DN")
DNname= data.frame(DNtype, p_DNid, DNpts)
names(DNname)[1] = "type"
names(DNname)[2] = "post_id"
names(DNname)[3] = "xyz"
# add specific type names if available
gs4_auth("ts587@cam.ac.uk")
DN_types = read_sheet("1M-2U0zV_T6BkdQZuZ2qVBGz6F_PD71MH0p5CJ71G38g", sheet = "matching")
DN_types_idsR = tibble(DN_types$DN_name, DN_types$RHS_newest_ID)
DN_types_idsL = tibble(DN_types$DN_name, DN_types$LHS_newest_ID)
names(DN_types_idsR)[1] = "type"
names(DN_types_idsR)[2] = "post_id"
names(DN_types_idsL)[1] = "type"
names(DN_types_idsL)[2] = "post_id"
DN_types_ids = bind_rows(DN_types_idsR,DN_types_idsL)
# i = 1
# while (i < nrow(DN_types_ids)) {
#     if (!flywire_islatest(DN_types_ids$post_id[i]))
#   {DN_types_ids$post_id[i] = flywire_latestid(DN_types_ids$post_id[i])}
#   else print(i)
#   i = i+1
# }
DN_types_ids = subset(DN_types_ids,(DN_types_ids$type != "-"))
DN_types_ids = subset(DN_types_ids,(DN_types_ids$post_id != "not found"))
DN_types_ids$type = as.character(DN_types_ids$type) 
DN_types_ids$post_id = as.character(DN_types_ids$post_id) 
DNname$type[match(DN_types_ids$post_id, DNname$post_id)] = DN_types_ids$type
# DN_Namiki = read_sheet("1vwlsODui1_kHCvLSyEWkbOoTwBZa52nJQNV_7BVSXc4", sheet = "ID_update")
# DN_namiki_ids = tibble(DN_Namiki$`neuron number`, DN_Namiki$flywire.id)
# DN_namiki_ids = subset(DN_namiki_ids,(!is.na(DN_namiki_ids$`DN_Namiki$flywire.id`)))
# names(DN_namiki_ids)[1] = "type"
# names(DN_namiki_ids)[2] = "post_id"
# DNname$type[match(DN_namiki_ids$post_id, DNname$post_id)] = DN_namiki_ids$type



# ANs
gs4_auth("tomke.stuerner@googlemail.com")
AN = read_sheet("10T0JE6nVSz_uUdoHGOpV2odO_k75-arRPKdOiBXlS80", sheet = "ANs")
ANpts= subset(AN$flywire.xyz,(!is.na(AN$flywire.xyz)))
ANid = flywire_xyz2id(ANpts, rawcoords=TRUE)
ANtype = replicate(length(ANid), "AN")
ANname= data.frame(ANtype, ANid, ANpts)
names(ANname)[1] = "type"
names(ANname)[2] = "post_id"
names(ANname)[3] = "xyz"
# matching sheet, 34448 matches from FAFB matching sheet
matches = read_sheet("1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw", sheet = "hemibrain")
matchespts= subset(matches$flywire.xyz,(!is.na(matches$flywire.xyz)))
matchesid = flywire_xyz2id(matchespts, rawcoords=TRUE)
matchestype = subset(matches$cell.type,(matches$flywire.xyz!= 'NA'))
matchesname= data.frame(matchestype, matchesid, matchespts)
matchesname= matchesname[complete.cases(matchesname), ]
names(matchesname)[1] = "type"
names(matchesname)[2] = "post_id"
names(matchesname)[3] = "xyz"

fmatches = read_sheet("1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw", sheet = "FAFB")
fmatchespts= subset(fmatches$flywire.xyz,(!is.na(fmatches$flywire.xyz)))
fmatchesid = flywire_xyz2id(fmatchespts, rawcoords=TRUE)
fmatchestype = subset(fmatches$cell.type,(fmatches$flywire.xyz!= 'NA'))
fmatchesname= data.frame(fmatchestype, fmatchesid, fmatchespts)
fmatchesname = subset(fmatchesname,(fmatchesname$fmatchestype!= 'NA'))
fmatchesname = subset(fmatchesname,(fmatchesname$fmatchestype!= "uncertain"))
fmatchesname= fmatchesname[complete.cases(fmatchesname), ]
names(fmatchesname)[1] = "type"
names(fmatchesname)[2] = "post_id"
names(fmatchesname)[3] = "xyz"
# Tomke neurons of interest
gs4_auth("ts587@cam.ac.uk")
NoI = read_sheet("1HcaY5dXD_SoUuQ2Y7dE9q1pc6XX6nsfLDl26a9MxCOE", sheet = "All_NoI")
NoIpts= subset(NoI$flywire.xyz,(!is.na(NoI$flywire.xyz)))
NoIid = flywire_xyz2id(NoIpts, rawcoords=TRUE)
NoItype = subset(NoI$name,(NoI$flywire.xyz!= 'NA'))
NoIname= data.frame(NoItype, NoIid, NoIpts)
NoIname= NoIname[complete.cases(NoIname), ]
names(NoIname)[1] = "type"
names(NoIname)[2] = "post_id"
names(NoIname)[3] = "xyz"
# put all together in one dataframe
flywire_types_all= bind_rows(MBONname, LHCENTname, DNname, ANname, matchesname, NoIname, fmatchesname)
# flywire_types= unique(flywire_types)
flywire_types= flywire_types_all%>% distinct(post_id, .keep_all = TRUE)
nrow(flywire_types)
saveRDS(flywire_types, file="flywire_types_05.10.Rds", version = NULL,
        compress = TRUE, refhook = NULL)

# 5239 typed flywire neurons
