#/bin/bash
OPTION=""
while [ -n "$1" ]
do
	case "$1" in
		--DGEMM_BLAS)
			opt="DGEMM_BLAS"
			OPTION="$OPTION $opt"
			;;
		--DGEMM_opt1)
			opt="DGEMM_opt1"
			OPTION="$OPTION $opt"
			;;
		--DGEMM_opt2)
			opt="DGEMM_opt2"
			OPTION="$OPTION $opt"
			;;
		--DGEMM_opt3)
			opt="DGEMM_opt3"
			OPTION="$OPTION $opt"
	esac
shift
done
while [ 1 ]
do
	echo "1. Запустить отлаженную программу"
	echo "2. Скомпилировать новую версию, сохраняя старую, или первая компиляция"
	echo "3. Перекомпилировать новую версию, без сохранения"
	echo "4. Выход"
	read choise
	
	case $choise in
		#Запуск отлаженой готовой программы
		1)
			#Проверка существования файла
			if [ -e $PWD/old_version.exe ]
			then
				./old_version.exe $OPTION
			else
				if [ -e $PWD/new_version.exe ]
				then	
					./new_version.exe $OPTION
				#Если ни один файл не найден
				else
					echo "File not found"
				fi
			fi
			;;
		#компиляция новой программы с сохранением старой
		2)
			echo "Enter name file.cpp"
			read file
			if [ -e $PWD/new_version.exe ]
			then
				#Если есть уже отлаженая программа
				#меняется её название, компилируется новая
				mv new_version.exe old_version.exe
				g++ -o new_version.exe $file
				./new_version.exe $OPTION
			else
				#если первая компиляция
				g++ -o new_version.exe $file
				./new_version.exe $OPTION
			fi
			;;
		3)
			#Перезапись предидущей компиляции
			echo "Enter name file.cpp"
			read file
			g++ -o new_version.exe $file
			./new_version.exe $OPTION
			;;
		4)
			exit 0
			;;
	esac
done
exit 0
