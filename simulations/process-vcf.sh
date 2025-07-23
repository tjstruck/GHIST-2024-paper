
declare -a arr=("GHIST-bottleneck" "GHIST-split-isolation" "GHIST-secondary-contact" "GHIST-admixture")
for i in "${arr[@]}"
do
    echo $i
    # Change sample names to something to make them more identifible
    bcftools reheader --samples  <(cut -f1 ${i}/${i}.popfile) ${i}/${i}.vcf | bcftools norm -d snps -o temp1-${i}.vcf
    # Make the reference allele the ancestral allele
    grep "^##" temp1-${i}.vcf > temp2-${i}.vcf
    echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele. Format: AA|REF|ALT|IndelType. AA: Ancestral allele, REF:Reference Allele, ALT:Alternate Allele, IndelType:Type of Indel (REF, ALT and IndelType are only defined for indels)">
' >> temp2-${i}.vcf
    grep "^#C" temp1-${i}.vcf >> temp2-${i}.vcf
    grep -v "^#" temp1-${i}.vcf | awk '{$8="AA="$4} 1' {O,}FS="\t" >> temp2-${i}.vcf
    mv temp2-${i}.vcf ${i}/${i}.vcf
    rm temp1-${i}.vcf
    gzip -9 ${i}/${i}.vcf
done


