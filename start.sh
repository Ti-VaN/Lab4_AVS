#/bin/bash
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
			echo "Enter parameters:"
			read param
			if [ -e $PWD/old_version.exe ]
			then
				./old_version.exe $param
			else
				if [ -e $PWD/new_version.exe ]
				then	
					./new_version.exe $param
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
			echo "Enter parameters:"
			read param
			if [ -e $PWD/new_version.exe ]
			then
				#Если есть уже отлаженая программа
				#меняется её название, компилируется новая
				mv new_version.exe old_version.exe
				g++ -o new_version.exe $file
				./new_version.exe $param
			else
				#если первая компиляция
				g++ -o new_version.exe $file
				./new_version.exe $param
			fi
			;;
		3)
			#Перезапись предидущей компиляции
			echo "Enter name file.cpp"
			read file
			echo "Enter parameters:"
			read param
			g++ -o new_version.exe $file
			./new_version.exe $param
			;;
		4)
			exit 0
			;;
	esac
done
exit 0
