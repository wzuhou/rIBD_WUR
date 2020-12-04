#! bin/bash
#----------------------------------------------------------------
#sh 0.nIBD_pair_pipeline_wind.sh DB GrMwB
#sh 0.nIBD_pair_pipeline_wind.sh GrMw GrMwB
#sh 0.nIBD_pair_pipeline_wind.sh DB FriFwB
#sh 0.nIBD_pair_pipeline_wind.sh FriFw FriFwB

for breed in `less list_group_cl.breed.large`; do \
sh 0.nIBD_pair_pipeline_wind.sh $breed  ${breed}B
sh 0.nIBD_pair_pipeline_wind.sh DB ${breed}B
sh 0.nIBD_pair_pipeline_wind.sh ${breed}B Jav
sh 0.nIBD_pair_pipeline_wind.sh ${breed}B  sebr
sh 0.nIBD_pair_pipeline_wind.sh ${breed}B  WelsumB
sh 0.nIBD_pair_pipeline_wind.sh ${breed}B jav_sebr.EIKENB
sh 0.nIBD_pair_pipeline_wind.sh ${breed}B jav_sebr
done

breed=KraiK
sh 0.nIBD_pair_pipeline_wind.sh $breed  ${breed}FwB
sh 0.nIBD_pair_pipeline_wind.sh DB ${breed}FwB
sh 0.nIBD_pair_pipeline_wind.sh ${breed}FwB Jav
sh 0.nIBD_pair_pipeline_wind.sh ${breed}FwB  sebr
sh 0.nIBD_pair_pipeline_wind.sh ${breed}FwB  WelsumB
sh 0.nIBD_pair_pipeline_wind.sh ${breed}FwB jav_sebr.EIKENB
sh 0.nIBD_pair_pipeline_wind.sh ${breed}FwB jav_sebr
#sh 0.nIBD_pair_pipeline_wind.sh DB_cl.neoBantam  jav_sebr.EIKENB
#sh 0.nIBD_pair_pipeline_wind.sh DB_cl.neoBantam sea

