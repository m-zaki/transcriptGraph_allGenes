# Loop function for the bedtools coverage

# 08/04/16

module load apps/bedtools/2.25.0/gcc-4.4.7 


source /home/zfadlullah/script/pbs/troodonHeader.sh

# Location of bed file (bed files which orginate from bam2bed)
bed_folder=/home/zfadlullah/tmp/GL26_NEB
bed_list=(`ls ${bed_folder}/*.bed`)

# Location of ensemble reference bed file
# Refer to link on how to generate the file -
EnsBED=/scratch/zaki/file/ref/zaki_customBED/EnsDb.Hsapiens.v75_proteinCoding_transcript.bed

# Bedtools Coverage output
perBase_folder=${bed_folder}/perBase
mkdir ${perBase_folder}


logDir=${perBase_folder}/log
mkdir ${logDir}

jname=bedCov_${startTime}

for file in ${bed_list[*]}
do
  f=`basename $file .bed`
  edited=${f}.txt
  
waitForQueue

  command="bedtools coverage -a ${EnsBED} -b ${file} -d -s > ${perBase_folder}/${edited}"

  createQsubScript "${command}" ${referrer}
  qsub -N ${jname} -e ${logDir}/${f}_bedCov.txt -o ${logDir}/bedCov_${f}_${startTime}_out.txt ${referrer}

done
waitForName ${jname}



source /home/zfadlullah/script/pbs/troodonFooter.sh

# -a = The EnsDb bed file generated from stage 1
# -b = the bed file which orginated from a BAM file
# -d = report per base coverage
##### Choose one of below depending on the experiment type
# -s = only report sequence with the same strand 
# -S = report opposite strandess (for duTP protocol)
