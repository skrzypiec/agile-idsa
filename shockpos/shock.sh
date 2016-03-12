echo "Path of data, i.e: /../data/data_s15/ "
read path

echo "Model name:"
read model

echo "Please enter bounce step:"
read bounce_step

echo "Please enter last step:"
read last_step

if [ $bounce_step -lt 10 ] && [ $last_step -gt 9 ]
then
	let n=9
	for ((i=bounce_step; i<=n; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_000$i.stp > temp_000$i.txt
	sed '1,11d' temp_000$i.txt > tempa_000$i.txt
	sed '104,209d' tempa_000$i.txt > temp_000$i.txt
	sed '5,220d' .${path}${model}_000$i.stp > time_000$i.txt
	sed '1,3d' time_000$i.txt > timea_000$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_000$i.txt > time_000$i.txt
	done
else
	if [ $bounce_step -lt 10 ] && [ $last_step -lt 9 ]
	then
	for ((i=bounce_step; i<=last_step; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_000$i.stp > temp_000$i.txt
	sed '1,11d' temp_000$i.txt > tempa_000$i.txt
	sed '104,209d' tempa_000$i.txt > temp_000$i.txt
	sed '5,220d' .${path}${model}_000$i.stp > time_000$i.txt
	sed '1,3d' time_000$i.txt > timea_000$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_000$i.txt > time_000$i.txt
	done
	fi
fi

if [ $bounce_step -lt 100 ] && [ $bounce_step -gt 9 ] && [ $last_step -gt 99 ]
then
	let n=99
	for ((i=bounce_step; i<=n; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_00$i.stp > temp_00$i.txt
	sed '1,11d' temp_00$i.txt > tempa_00$i.txt
	sed '104,209d' tempa_00$i.txt > temp_00$i.txt
	sed '5,220d' .${path}${model}_00$i.stp > time_00$i.txt
	sed '1,3d' time_00$i.txt > timea_00$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_00$i.txt > time_00$i.txt
	done
else
	if [ $bounce_step -lt 100 ] && [ $bounce_step -gt 9 ] && [ $last_step -lt 99 ]
	then
	for ((i=bounce_step; i<=last_step; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_00$i.stp > temp_00$i.txt
	sed '1,11d' temp_00$i.txt > tempa_00$i.txt
	sed '104,209d' tempa_00$i.txt > temp_00$i.txt
	sed '5,220d' .${path}${model}_00$i.stp > time_00$i.txt
	sed '1,3d' time_00$i.txt > timea_00$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_00$i.txt > time_00$i.txt
	done
	fi
fi


if [ $bounce_step -lt 9 ] && [ $last_step -gt 99 ]
then
	let n=99
	for ((i=10; i<=n; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_00$i.stp > temp_00$i.txt
	sed '1,11d' temp_00$i.txt > tempa_00$i.txt
	sed '104,209d' tempa_00$i.txt > temp_00$i.txt
	sed '5,220d' .${path}${model}_00$i.stp > time_00$i.txt
	sed '1,3d' time_00$i.txt > timea_00$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_00$i.txt > time_00$i.txt
	done
else
	if [ $bounce_step -lt 9 ] && [ $last_step -lt 99 ]
	then
	for ((i=10; i<=last_step; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_00$i.stp > temp_00$i.txt
	sed '1,11d' temp_00$i.txt > tempa_00$i.txt
	sed '104,209d' tempa_00$i.txt > temp_00$i.txt
	sed '5,220d' .${path}${model}_00$i.stp > time_00$i.txt
	sed '1,3d' time_00$i.txt > timea_00$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_00$i.txt > time_00$i.txt
	done
	fi
fi


if [ $bounce_step -lt 1000 ] && [ $bounce_step -gt 99 ] && [ $last_step -lt 1000 ]
then
	for ((i=bounce_step; i<=last_step; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_0$i.stp > temp_0$i.txt
	sed '1,11d' temp_0$i.txt > tempa_0$i.txt
	sed '104,209d' tempa_0$i.txt > temp_0$i.txt
	sed '5,220d' .${path}${model}_0$i.stp > time_0$i.txt
	sed '1,3d' time_0$i.txt > timea_0$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_0$i.txt > time_0$i.txt
	done
else
	if [ $bounce_step -lt 100 ] && [ $last_step -lt 1000 ]
	then
	for ((i=100; i<=last_step; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_0$i.stp > temp_0$i.txt
	sed '1,11d' temp_0$i.txt > tempa_0$i.txt
	sed '104,209d' tempa_0$i.txt > temp_0$i.txt
	sed '5,220d' .${path}${model}_0$i.stp > time_0$i.txt
	sed '1,3d' time_0$i.txt > timea_0$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_0$i.txt > time_0$i.txt
	done
	fi
fi

if [ $bounce_step -lt 1000 ] && [ $bounce_step -gt 99 ] && [ $last_step -gt 999 ]
then
	let n=999
	for ((i=bounce_step; i<=n; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_0$i.stp > temp_0$i.txt
	sed '1,11d' temp_0$i.txt > tempa_0$i.txt
	sed '104,209d' tempa_0$i.txt > temp_0$i.txt
	sed '5,220d' .${path}${model}_0$i.stp > time_0$i.txt
	sed '1,3d' time_0$i.txt > timea_0$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_0$i.txt > time_0$i.txt
	done
else
	if  [ $bounce_step -lt 100 ] && [ $last_step -gt 999 ]
	then
	let n=999
	for ((i=100; i<=n; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_0$i.stp > temp_0$i.txt
	sed '1,11d' temp_0$i.txt > tempa_0$i.txt
	sed '104,209d' tempa_0$i.txt > temp_0$i.txt
	sed '5,220d' .${path}${model}_0$i.stp > time_0$i.txt
	sed '1,3d' time_0$i.txt > timea_0$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_0$i.txt > time_0$i.txt
	done
	fi
fi

if [ $bounce_step -gt 999 ]
then
	for ((i=bounce_step; i<=last_step; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_$i.stp > temp_$i.txt
	sed '1,11d' temp_$i.txt > tempa_$i.txt
	sed '104,209d' tempa_$i.txt > temp_$i.txt
	sed '5,220d' .${path}${model}_$i.stp > time_$i.txt
	sed '1,3d' time_$i.txt > timea_$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_$i.txt > time_$i.txt
	done
else
	if [ $bounce_step -lt 1000 ] && [ $last_step -gt 999 ]
	then
	for ((i=1000; i<=last_step; i++))
	do
	awk '{$1=""; $2=""; $5=""; $6=""; $7=""; $8=""; $9=""; $10=""; $11=""; $12=""; $13=""; sub("  ", " "); print}' .${path}${model}_$i.stp > temp_$i.txt
	sed '1,11d' temp_$i.txt > tempa_$i.txt
	sed '104,209d' tempa_$i.txt > temp_$i.txt
	sed '5,220d' .${path}${model}_$i.stp > time_$i.txt
	sed '1,3d' time_$i.txt > timea_$i.txt
	awk '{$1=""; sub("  ", " "); print}' timea_$i.txt > time_$i.txt
	done
	fi
fi





rm tempa*

######################## wycinanie dwóch kolumn z plików w zależności od bounce i last step

rm timea_*
cat time_* > times_1.txt

let j=$last_step-$bounce_step

if [ $bounce_step -lt 10 ]
then

	for ((i=0; i<=j; i++))
	do
	cat time_000$bounce_step.txt > timea_$i.txt
	done
fi

if [ $bounce_step -lt 100 ] && [ $bounce_step -gt 9 ]
then
	for ((i=0; i<=j; i++))
	do
		cat time_00$bounce_step.txt > timea_$i.txt
	done
fi

if [ $bounce_step -lt 1000 ] && [ $bounce_step -gt 99 ]
then
	for ((i=0; i<=j; i++))
	do
		cat time_0$bounce_step.txt > timea_$i.txt
	done
fi

if [ $bounce_step -gt 999 ]
then
	for ((i=0; i<=j; i++))
	do
		cat time_$bounce_step.txt > timea_$i.txt
	done
fi


cat timea_* > times_2.txt

rm time_*
rm timea_*

paste times_1.txt times_2.txt | awk '{ print $1 - $2 }' > time.txt

rm times_*


################################# time issues

if [ $bounce_step -gt 10 ]
then
    let f=10
else
	let f=$bounce_step
fi

if [ $bounce_step -gt 100 ] 
then
    let g=100
	elif [ $bounce_step -gt 9 ]
	then
		let g=$bounce_step
	else 
		let g=10
fi


if [ $bounce_step -gt 1000 ]
then
    let h=1000
	elif [ $bounce_step -gt 99 ]
	then
		let h=$bounce_step
	else 
		let h=100
fi


if [ $bounce_step -gt 9999 ]
then
    let k=9999
	elif [ $bounce_step -gt 999 ]
	then
		let k=$bounce_step
	else 
		let k=1000
fi


if [ $last_step -gt 9 ]
then
    let q=9
else
	let q=$last_step
fi

if [ $last_step -gt 99 ]
then
    let p=99
else
	let p=$last_step
fi


if [ $last_step -gt 999 ]
then
    let r=999
else
	let r=$last_step
fi

if [ $last_step -gt 9999 ]
then
    let s=9999
else
	let s=$last_step
fi

if [ $bounce_step -lt 10 ] 
then

	for ((i=f; i<=q; i++))
	do
	./finddiv temp_000${i}.txt; 
	done
fi

if [ $bounce_step -lt 100 ] 
then
	for ((i=g; i<=p; i++))
	do
		./finddiv temp_00${i}.txt; 
	done
fi

if [ $bounce_step -lt 1000 ] 
then
	for ((i=h; i<=r; i++))
	do
		./finddiv temp_0${i}.txt; 
	done
fi

	for ((i=k; i<=s; i++))
	do
		./finddiv temp_${i}.txt; 
	done


cat temp_* > rad.txt

#rm temp_*

pr -mts"    " time.txt rad.txt > shockpos.txt

rm time.txt
rm rad.txt

output=${model}_${bounce_step}.png

gnuplot -e "outputname='${output}'" 'wykres.plt' 





