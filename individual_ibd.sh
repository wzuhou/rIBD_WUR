#! /bin/bash
for id in `less ./List_animal_ID/List.GrMwB`;\
do grep $id ibdpair.GrMwB_DB_sort >GrMwB_DB.${id};\
grep $id ibdpair.GrMw_GrMwB_sort> GrMw_GrMwB.${id};\
awk '$6>=33080000 && $6<35370000' GrMwB_DB.${id}>GrMwB_DB.${id}.FK1M;\
awk '$6>=33080000 && $6<35370000' GrMw_GrMwB.${id}>GrMw_GrMwB.${id}.FK1M;\
done

