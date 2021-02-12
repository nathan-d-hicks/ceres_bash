while read p;
do
	sbatch ~/bin/download_sra.sh $p
	sleep 1
done < $1
