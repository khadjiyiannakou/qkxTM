while read -r LINE
do
    list=( $LINE )
    numElements=${#list[@]}

numElements=`echo "$numElements/2" | bc `
    
    echo "
{
L=$numElements
for(i = 0 ; i < L ; i++)
{
if(\$((i+1)*2) <= 0){printf \"%+16.15e - %16.15ei\t\",\$((i*2)+1),sqrt((\$((i+1)*2))*(\$((i+1)*2)))}
if(\$((i+1)*2) > 0){printf \"%+16.15e + %16.15ei\t\",\$((i*2)+1),sqrt((\$((i+1)*2))*(\$((i+1)*2)))}
}
printf \"\n\"
}
" > awk.txt
echo $LINE | awk -f awk.txt

done < 2.dat