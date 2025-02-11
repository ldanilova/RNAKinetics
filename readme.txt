Ключи

"-id" = "--id" - идентификатор задачи
"-seq" = "--seq" - последовательность (все равно какую, ДНК или РНК, если ДНК, то я все t меняю на u, и дальше работаю как с РНК). Минимальную длину последовательности можно сделать 20.
"-tm" = "--time" - время (по умолчанию 3)
"-gr" = "--growth" - константа скорости роста (по умолчанию 40)
"-ncl" = "--nucleation" - константа скорости образования очередной пары (по уолчанию 10^7)
"-m" = "--m" - количество проходов (по умолчанию 10)
"-in" = "--infile" - входной файл с последовательностью в fasta-формате (в качестве имени можно дать id). В принципе, я умею обрабатывать несколько последовательностей по очереди, но сейчас будет браться только первая. Потом можно будет переделать, если нужно.
"-out" = "--outfile" - выходной файл. Если его не будет, то результаты запишутся в файл "id".out
"-hl" = "--helixLength" - минимальная длина спирали (по умолчанию 3)
"-gu" = "--gu" - учитываем или нет концевые GU-пары (по умолчанию false)
"-ht" = "--helixThreshold" - порог отсечения спиралей по энергии спирали (по умолчанию -7000)
"-td" = "--threshold_dis" - порог на константу разрушения при объединении в группу (по умолчанию 1.e+5)
"-tf" = "--threshold_form" - порог на константу образования при объединении в группу (по умолчанию 1.e+5)



Java-ключи

-Xmxn Specify the maximum size, in bytes, of the memory allocation pool. This value must a multiple of 1024 greater than 2MB. Append the letter k or K to indicate kilobytes, or m or M to indicate megabytes. The default value is 64MB. The upper limit for this value will be approximately 4000m on Solaris 7 and Solaris 8 SPARC platforms and 2000m on Solaris 2.6 and x86 platforms, minus overhead amounts. Examples:        
       -Xmx83886080
       -Xmx81920k
       -Xmx80m



=========================================
Запуск из командной строки

из папки classes

java rna1.rna1 -id 112 -in test.txt

java rna1.rna1 -id 112 -seq GCCGCCGUAGCUCAGCCCGGGAGAGCGCCCGGCUGAAGACCGGGUUGUCCGGGGUUCAAGUCCCCGCGGCGGCA

запуск jar-архива

java -jar rna1.jar -id 111 -in test.txt

java -jar rna1.jar -id 111 -seq GCCGCCGUAGCUCAGCCCGGGAGAGCGCCCGGCUGAAGACCGGGUUGUCCGGGGUUCAAGUCCCCGCGGCGGCA

==========================================
Выходной файл в XML формате
<RNAKinetics> - главный тег с атрибутом ID
-----------------------------
Input_Params -  входные параметры программы

Sequence - последовательность с длиной length 
Time - время (t)
M - количество проходов
k_growth- константа роста цепи
k_nucl - константа образования оч. пары
helix_threshold - порог отсечения спиралей
destroy_threshold - порог на константу разрушения при сборке группы
formation_threshold - порого на константу образования
min_helix_length - минимальная длина спиралей
allows_final_GU - учитываем или нет концевые GU-пары
structures_graph_step - количество шагов для графиков структур и спиралей (graph step for structures). Сейчас это 100.

<Input_Params>
	<Sequence length = "76">ggggcuauagcucagcugggagagcgccugcuuugcacgcaggaggucugcgguucgaucccgcauagcuccacca</Sequence>
	<Time>5.0</Time>
	<M>100</M>
...
</Input_Params>


------------------------------------------
Statistics - статистика работы программы

permutations - количество структурных перестановок за все время работы программы
grow_time - момент времени, когда последовательность выросла (grown_time) или если она не успела вырасти, то позицию, до которой успела (grown_position).
work_time - время работы программы

<Statistics>
	<permutations>23419</permutations>
	<grow_time>2.385229690183856</grow_time>
	<work_time>1062</work_time>
...
</Statistics>


---------------------------------
Heliset - Список спиралей состоит из Helix 

Helix - спираль
атрибуты:
ID - номер в Heliset
Energy - энергия спирали
Length - длина спирали
BreakConst - константа скорости разрушения спирали

<Graph> - график зависимости вероятности реализации спирали от времени. Состоит из point
LeftShoulder - левое плечо спирали с координатами начала и конца
RightShoulder - перевернутое правое плечо с координатами


<Helix ID = "0" Energy = "-7600.0" Length = "4" BreakConst = "131.42609948192765">
		<Graph>
			<point x = "0" >0.0</point>
			<point x = "1" >0.0</point>
		 	<point x = "2" >0.040791579383376436</point>
		 	<point x = "3" >0.11363065935476988</point>
		 	<point x = "4" >6.99203347425037E-4</point>
		 	<point x = "5" >0.0</point>
		 	<point x = "6" >0.0</point>
			...
		 </Graph>
		<LeftShoulder start = "2" end = "5" >ggcu</LeftShoulder>
		<RightShoulder start = "16" end = "13" >ucga</RightShoulder >		
	</Helix>


------------------------------------
GroupList - список групп, состоит из Group 

Size - размер группы
Group - группа спиралей состоит из GroupElement 
Атрибуты: 
ID - номер в списке групп
Completed - полная группа или нет
BestElementID - номер лучшего элемента

GroupElement - элементы группы, т.е. какая-то структура
Атрибуты:
ID - номер в группе

Energy - энергия структуры
LoopsEnergy - энергия петель структуры
Folding - фолдинг структуры в точечно-скобочной записи
SlidePairs - скользящие пары, состоят из перечечения Cross и альтернативных спариваний Alter
HelixID - номера спиралей из списка спиралей, которые образуют этот фолдинг

<GroupList>
	<Size>260</Size>
	<Group ID = "216" Completed = "true" BestElementID = "1" >
		<GroupElement ID = "0" IsBest = "false" >
			<Energy>-17400.0</Energy>
			<LoopsEnergy>20200.0</LoopsEnergy>
			<Folding>.............(((..(((.....................................)))......)))......</Folding>
			<SlidePairs>
				<Cross>16</Cross>
				<Alter>66</Alter>
				<Alter>62</Alter>
			</SlidePairs>
			<SlidePairs>
				<Cross>17</Cross>
				<Alter>65</Alter>
				<Alter>61</Alter>
			</SlidePairs>
			<HelixID>10</HelixID>
			<HelixID>13</HelixID>			
		</GroupElement>
	</Group> 
</GroupList>

--------------------------------------------------
<Results> -  результаты работы программы состоит из ResultElement

ResultElement - результат с атрибутами Time и GroupID - номером группы из списка групп
Из этой группы надо взять элемент с атрибутом IsBest = "true"

Graph - график вероятности реализации структуры от времени. Состоит из point

<Results>
	<ResultElement Time = "3.3378043" GroupID = "5" >
		 <Graph>
			<point x = "0" >0.0</point>
			<point x = "1" >0.0</point>
		 	<point x = "2" >0.040791579383376436</point>
		 	<point x = "3" >0.11363065935476988</point>
		 	<point x = "4" >6.99203347425037E-4</point>
		 	<point x = "5" >0.0</point>
		 	<point x = "6" >0.0</point>
			...
		 </Graph>

	</ResultElement>
</Results>



=============================================
Если в последовательности нет спиралей, то пишется 
<Message>There is no helices</Message>

--------------------------------------------

Если в последовательности только одна спираль, то пишется:
<Message>There is only one helix</Message> , 
Heliset и 

==============================================
log-файл

записывается текущий запуск (run) программы, общее число шагов и текущее время в секундах 

current run m = 10
M = 100
start time = 104458511
current time = 110795104

Если не хватило памяти на каком-то запуске, то еще припишется "Out of memory" и программа закончит свою работу.

Log-файл обновляется когда начинается новый запуск. 

После формирования output-файла с результатами в log-файл записывается значение ключа "-m".



Так же в log-файл может быть записано следующее:
"There is no helices" -  если в последовательнсти нет ниодной спирали
"There is only one helix" - если в последовательности только одна спираль
"Number of next step doesn't find" -  если все нормально, то этой записи быть не должно. Она означает, что по каким-то причинам (даже не могу придумать по каким, но на всякий случай сделала) мы не можем определить что должно произойти за следующий шаг: вырасти последовательность или разрушиться спираль или образоваться спираль.
"File is not exist" - при попытке открыть несуществующий файл
"No data!"- если в открываемом файле нет данных
"Heliset creation. Out of memory" - если не хватило памяти при создании списка потенциальных спиралей
"Out of memory" - если не хватило памяти на других этапах работы программы
"Uncaught Exception" - если произошло еще что-то.