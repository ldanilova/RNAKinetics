Keys

"-seq" = "--seq" - DNA or RNA sequence. The length of sequence should be more than 20 
"-in" = "--infile" - Input file with sequence(s) (<seq ID> tab <sequence>).
"-id" = "--id" - ID of task ("test" by default)
"-tm" = "--time" - Final time (3 by default)
"-gr" = "--growth" - Chain's growth rate (40 by default)
"-ncl" = "--nucleation" - The kinetic constant of one marginal complementary pair locking (Nucleation constant) (1e+7 by default)
"-m" = "--m" - Number of runs (100 by default)
"-out" = "--outfile" - Output file (<id>.xml by default)
"-hl" = "--helixLength" - Min helix length (4 by default)
"-gu" = "--gu" - Allow terminal GU. 0 is false and 1 is true. False by default.
"-ht" = "--helixThreshold" - Helix energy threshold (-7000 by default)
"-td" = "--threshold_dis" - Decay kinetic constant for grouping (1.e+5 by default)
"-tf" = "--threshold_form" - Formation kinetic constant for grouping (1.e+5 by default)

Java-keys for increasing memory

-Xmxn Specify the maximum size, in bytes, of the memory allocation pool. This value must a multiple of 1024 greater than 2MB. Append the letter k or K to indicate kilobytes, or m or M to indicate megabytes. The default value is 64MB. The upper limit for this value will be approximately 4000m on Solaris 7 and Solaris 8 SPARC platforms and 2000m on Solaris 2.6 and x86 platforms, minus overhead amounts. Examples:        
       -Xmx83886080
       -Xmx81920k
       -Xmx80m

=========================================
Example of command line 

java -Xmx200M -jar rna6.jar -in test.txt

java -jar rna6.jar -id 111 -seq GCCGCCGUAGCUCAGCCCGGGAGAGCGCCCGGCUGAAGACCGGGUUGUCCGGGGUUCAAGUCCCCGCGGCGGCA

==========================================
Output XML file tags:

RNAKinetics - main tag with ID attribute
---------------------------------------
Input_Params tag -  Input parameters:
  Sequence tag - input sequence with 'length' attribute 
  Time tag - final time
  M tag - Number of runs
  k_growth tag - Chain's growth rate
  k_nucl tag - Nucleation constant
  helix_threshold tag - Helix energy threshold
  destroy_threshold tag - Decay kinetic constant for grouping
  formation_threshold tag - Formation kinetic constant for grouping
  min_helix_length tag - Min helix length
  allows_final_GU tag - Allow terminal GU
  structures_graph_step tag - Graph step for structures

<Input_Params>
	<Sequence length = "76">ggggcuauagcucagcugggagagcgccugcuuugcacgcaggaggucugcgguucgaucccgcauagcuccacca</Sequence>
	<Time>5.0</Time>
	<M>100</M>
...
</Input_Params>


------------------------------------------
Statistics tag - Statistics of program run

  permutations tag - Number of permutation during program run
  grow_time tag - Timestamp when the sequence has grown
  grown_position tag - Position in the sequence where it grew if the sequence hadn't grown by the final time
  work_time tag - Computer time of program run

<Statistics>
	<permutations>23419</permutations>
	<grow_time>2.385229690183856</grow_time>
	<work_time>1062</work_time>
...
</Statistics>

---------------------------------
Heliset -List of helices 

  Helix tag:
    atributes:
    ID - index in Heliset
    Energy - helix energy
    Length - helix length
    BreakConst - Decay kinetic constant 

  LeftShoulder tag - left shoulder of helix with the start and end coordinats
  RightShoulder tag - right shoulder of helix with the start and end coordinats
  Graph tag - plot of probability distribution of helix appearance as a function of Ò    

<Helix ID = "0" Energy = "-7600.0" Length = "4" BreakConst = "131.42609948192765">
		<Graph>
			0E0, 0E0, 0E0, 0E0, 0E0, 0E0, 0E0, 0E0, 8.14328E-4, 6.58596E-3, 1.01853E-2, 1.61663E-2,
		 </Graph>
		<LeftShoulder start = "2" end = "5" >ggcu</LeftShoulder>
		<RightShoulder start = "16" end = "13" >ucga</RightShoulder >		
</Helix>

------------------------------------
GroupList tag - list of groups which contains Group
 Size tag - group size
 Group tag - group which contains GroupElement. The first element in the group is the best.
  attribute: 
  ID - index in the group list
  Completed - is it completed group or not
  
 GroupElement tag - element of group is a structure.
  attribute:
  ID - index in group

  Energy tag - structure energy
  Folding tag - structure in dot bracket format
  SlidePairs tag - sliding pairs which contain intersection in the Cross tag and alternative pairs in the Alter tag
  HelixID tag - index of helix from list of helices which forms this structure

<GroupList>
	<Size>260</Size>
	<Group ID = "216" Completed = "true">
		<GroupElement ID = "0" >
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
Results tag - results of program run which contains ResultElement tags

ResultElement tag:
  attributes:
  Time - lifetime
  GroupID - group index in the list of groups

Graph tag - plot of probability distribution of structure appearance as a function of Ò

<Results>
	<ResultElement Time = "3.3378043" GroupID = "5" >
		 <Graph>
			0E0, 0E0, 0E0, 0E0, 0E0, 0E0, 0E0, 0E
			...
		 </Graph>

	</ResultElement>
</Results>

=============================================
If the input sequence doesn't have a helix, the following message is written:
<Message>There is no helices</Message>

--------------------------------------------
If the input sequence has the only helix, the following message is written:
<Message>There is only one helix</Message> , 
and Heliset with one helix

==============================================
log-file

Store current run of the program, number of runs, start time and current time in seconds, and errors messages (if any)  

current run m = 10
M = 100
start time = 104458511
current time = 110795104
