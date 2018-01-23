#Grep out contained events
zcat reviseAnnotations/naive.permuted.txt.gz | grep contained | gzip > txrevise_contained/naive.permuted.txt.gz
zcat reviseAnnotations/IFNg.permuted.txt.gz | grep contained | gzip > txrevise_contained/IFNg.permuted.txt.gz
zcat reviseAnnotations/SL1344.permuted.txt.gz | grep contained | gzip > txrevise_contained/SL1344.permuted.txt.gz
zcat reviseAnnotations/IFNg_SL1344.permuted.txt.gz | grep contained | gzip > txrevise_contained/IFNg_SL1344.permuted.txt.gz

#Copy full summary stats
cp -r reviseAnnotations/sorted/ txrevise_contained/ &

##For acLDL data
#Grep out contained events
zcat reviseAnnotations/Ctrl.permuted.txt.gz | grep contained | gzip > txrevise_contained/Ctrl.permuted.txt.gz
zcat reviseAnnotations/AcLDL.permuted.txt.gz | grep contained | gzip > txrevise_contained/AcLDL.permuted.txt.gz

#Copy full summary stats
cp -r reviseAnnotations/sorted/ txrevise_contained/ &
