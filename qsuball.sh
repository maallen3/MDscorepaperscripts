
#id=/scratch/Shares/dowell/counts/human20151211-182056/
id=/scratch/Shares/dowell/counts/human20160929-114122/

for pathandfilename in `ls $id*SRR2961*.sh`; do
echo $pathandfilename
qsub $pathandfilename
sleep 1
done
