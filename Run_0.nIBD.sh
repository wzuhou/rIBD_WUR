#! bin/bash
sh 0.nIBD_pair_pipeline_wind.sh DB DB_cl.neoBantam
sh 0.nIBD_pair_pipeline_wind.sh DB_cl.large DB_cl.neoBantam


sh 0.nIBD_pair_pipeline_wind.sh jav_sebr jav_sebr_cl.neoBantam
sh 0.nIBD_pair_pipeline_wind.sh jav_sebr_cl.large jav_sebr_cl.neoBantam


sh 0.nIBD_pair_pipeline_wind.sh jav_sebr.EIKENB jav_sebr_cl.neoBantam.EIKENB
sh 0.nIBD_pair_pipeline_wind.sh jav_sebr_cl.large jav_sebr_cl.neoBantam.EIKENB

sh 0.nIBD_pair_pipeline_wind.sh sea sea_cl.neoBantam
sh 0.nIBD_pair_pipeline_wind.sh sea_cl.large sea_cl.neoBantam

sh 0.nIBD_pair_pipeline_wind.sh full_bantam full_neoBantam
sh 0.nIBD_pair_pipeline_wind.sh full_large full_neoBantam

sh 0.nIBD_pair_pipeline_wind.sh full_bantam.EIKENB full_neoBantam.EIKENB
sh 0.nIBD_pair_pipeline_wind.sh full_large full_neoBantam.EIKENB
#----------------------------------------------------------------
sh 0.nIBD_pair_pipeline_wind.sh DB GrMwB
sh 0.nIBD_pair_pipeline_wind.sh GrMw GrMwB
sh 0.nIBD_pair_pipeline_wind.sh DB FriFwB
sh 0.nIBD_pair_pipeline_wind.sh FriFw FriFwB

sh 0.nIBD_pair_pipeline_wind.sh DB_cl.neoBantam  jav_sebr.EIKENB
sh 0.nIBD_pair_pipeline_wind.sh DB_cl.neoBantam sea

