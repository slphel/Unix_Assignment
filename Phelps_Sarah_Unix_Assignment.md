#UNIX Assignment Sarah Phelps

##Data Inspection

###Attributes of `fang_et_al_genotypes`

```
wc fang_et_al_genotypes.txt           # number of lines, words, and bytes
#   2783  2744038 11051939 fang_et_al_genotypes.txt 

du -h fang_et_al_genotypes.txt         # how large is the file
# 6.7M	fang_et_al_genotypes.txt

 awk -F "\t" '{print NF; exit}' fang_et_al_genotypes.txt        # number of columns
 # 986
 
```

By inspecting this file I learned that:

1. The file size is 6.7M
2. There are 2783 lines, 2744038 words, and 11051939 bytes in the data set
3. There are 986 columns in the file. This explains why the "head" and "tail"
commands are not the most useful.



### Attributes of `snp_position.txt`

```
du -h snp_position.txt 			 # inspecting file size
# 49K	snp_position.txt

awk -F "\t" '{print NF; exit}' snp_position.txt 		# how many columns are in the file
# 15

wc snp_position.txt 
# 984 13198 82763 snp_position.txt  # number of lines, words, and bytes

head -n 1 snp_position.txt 				# Read the header of the file

cat snp_position.txt			# Interactive view of the file


```

By inspecting this file I learned that:

1. The file is 49k
2. There are 15 columns in the file. The title of the columns are: SNP_ID, cdv_marker_id, Chromosome, 
Position, alt_pos, mult_positions, amplicon, cdv_map_feature.name, gene, candidate/random
Genaissance_daa_id, Sequenom_daa_id, count_amplicons, count_cmf, count_gene.
3. There are 984 lines in the data file. Since there is a header, we know there is 983 lines of data.



##Data Processing

###Maize Data

```
# Cleaning, Sorting, and Transposing Maize groups

awk '$3 ~ /ZMMIL|ZMMLR|ZMMMR/' fang_et_al_genotypes.txt > hi.txt	## Filter 3 maize groups to a new file

head -n 1 fang_et_al_genotypes.txt > header.txt 	## Add the header from the original fang file to a new file

cat header.txt hi.txt > combined-maize  	## add the filtered file containing the three maize groups and the file with the header to a new file

head -n 3 combined-maize ## check file

awk -F "\t" '{print NF; exit}' combined-maize ## print the number of colums in the file

cut -f 3-986 combined-maize > clean-maize 	## clean the file up with only necessary columns 

awk -f transpose.awk clean-maize > transposed-maize.txt 	## transpose cleaned file

sort -k1,1 transposed-maize.txt > sort-trans-maize 			## sort the transposed file by chromosome.


 
 
```

Here is my brief description of what this code does

This code started by filtering the three SNP groups of interest from the Fang_et_al_genotypes file.  
Then I had to add the header from the original fang file. Then I cut the file to only the three columns I needed. 
Using the code provided, the file was transposed and then sorted by ascending order of chromosome position.




###Teosinte Data

Cleaning, Sorting, and Transposing Teosinte groups

```

awk '$3 ~ /ZMPBA|ZMPIL|ZMPJA/' fang_et_al_genotypes.txt > bye.txt	 ## Filter 3 maize groups to a new file

cat header.txt bye.txt > combined-teos.txt 		## Add the header from the original fang file to a new file

head -n 3 combined-teos.txt		 ## check file

awk -F "\t" '{print NF; exit}' combined-teos.txt  ## print the number of colums in the file

cut -f 3-986 combined-teos.txt	> clean-teos.txt ## clean the file up with only necessary columns 

awk -f transpose.awk clean-teos.txt > transposed-teos.txt 		## transpose cleaned file

sort -k1,1 transposed-teos.txt > sort-trans-teos.txt			## sort the transposed file by chromosome.

```

Here is my brief description of what this code does

This code performed the same tasks on the teosinte groups as the previous code did on the maize groups. 



# Joining the Maize and teosinte file with the SNP file (2 separate files)
I included short descriptions on each line to state why I ran each code.
```

mv sort-trans-maize sort-trans-maize.txt					 ## rename as a txt file 

awk '{print $1 "\t" $3 "\t" $4}' snp_position.txt > cut.snp.txt			 ## cut the snp file to the 3 columns I need

head -n 1 cut.snp.txt 		#check file 

tail -n +2 cut.snp.txt | sort -k1,1 > sorted_cut.snp.txt 		## sort snp file to join

join -1 1 -2 1 sorted_cut.snp.txt sort-trans-maize.txt > join1.txt 		## join maize file to snp file 

join -1 1 -2 1 sorted_cut.snp.txt sort-trans-teos.txt > joint.txt 		## join teos file to snp file 



```

Here is my brief description of what this code does
This code cut the three columns I need from the snp file and created a new file. That new file was joined with the
transposed maize file to create the final file to manipulate the maize data. The cut SNP file was also joined with 
the teosinte file to create final file to manipulate the teosinte data.



# Creating new files for the maize data
## Descriptions are placed on the code line as needed. 

```
sort -k3,3 join1.txt > sjoin1.txt 		## sort by colum 3 (position)

awk '$2 == "1" { print $0 }' sjoin1.txt > chr1.maize.txt 		##chr1

awk '$2 == "2" { print $0 }' sjoin1.txt > chr2.maize.txt 		##chr2

awk '$2 == "3" { print $0 }' sjoin1.txt > chr3.maize.txt

awk '$2 == "4" { print $0 }' sjoin1.txt > chr4.maize.txt

awk '$2 == "5" { print $0 }' sjoin1.txt > chr5.maize.txt

awk '$2 == "6" { print $0 }' sjoin1.txt > chr6.maize.txt

awk '$2 == "7" { print $0 }' sjoin1.txt > chr7.maize.txt

awk '$2 == "8" { print $0 }' sjoin1.txt > chr8.maize.txt

awk '$2 == "9" { print $0 }' sjoin1.txt > chr9.maize.txt

awk '$2 == "10" { print $0 }' sjoin1.txt > chr10.maize.txt

sed 's/?/-/g' sjoin1.txt > sjoin1d.txt 	##replace missing values "?" with "-"



## Decreasing Position 

sort -k3,3nr join1.txt >  dm.txt			##Sort by decreasing position

awk '$2 == "1" { print $0 }' dm.txt > dchr1.maize.txt 		##chr1

awk '$2 == "2" { print $0 }' dm.txt > dchr2.maize.txt 		##chr2

awk '$2 == "3" { print $0 }' dm.txt > dchr3.maize.txt

awk '$2 == "4" { print $0 }' dm.txt > dchr4.maize.txt

awk '$2 == "5" { print $0 }' dm.txt > dchr5.maize.txt

awk '$2 == "6" { print $0 }' dm.txt > dchr6.maize.txt

awk '$2 == "7" { print $0 }' dm.txt > dchr7.maize.txt

awk '$2 == "8" { print $0 }' dm.txt > dchr8.maize.txt

awk '$2 == "9" { print $0 }' dm.txt > dchr9.maize.txt

awk '$2 == "10" { print $0 }' dm.txt > dchr10.maize.txt




## unknown positions 

awk '$3 == "unknown" { print $0 }' dm.txt > uk.maize.txt ##replace missing values "?" with "-"



## move all files to a maize folder 

mkdir maizefolder

mv chr1.maize.txt chr2.maize.txt chr3.maize.txt chr4.maize.txt chr5.maize.txt chr6.maize.txt chr7.maize.txt chr8.maize.txt chr9.maize.txt chr10.maize.txt maizefolder/

mv dchr1.maize.txt dchr2.maize.txt dchr3.maize.txt dchr4.maize.txt dchr5.maize.txt dchr6.maize.txt dchr7.maize.txt dchr8.maize.txt dchr9.maize.txt dchr10.maize.txt maizefolder/

mv uk.maize.txt maizefolder/

```



# Creating new files for the Teosinte Data.
##  Descriptions are placed on the code line as needed. 

```
sort -k3,3 joint.txt > sjoint.txt 	## sort by colum 3 (position)

awk '$2 == "1" { print $0 }' sjoint.txt > chr1.teos.txt 		##chr1

awk '$2 == "2" { print $0 }' sjoint.txt > chr2.teos.txt 		##chr2

awk '$2 == "3" { print $0 }' sjoint.txt > chr3.teos.txt

awk '$2 == "4" { print $0 }' sjoint.txt > chr4.teos.txt

awk '$2 == "5" { print $0 }' sjoint.txt > chr5.teos.txt

awk '$2 == "6" { print $0 }' sjoint.txt > chr6.teos.txt

awk '$2 == "7" { print $0 }' sjoint.txt > chr7.teos.txt

awk '$2 == "8" { print $0 }' sjoint.txt > chr8.teos.txt

awk '$2 == "9" { print $0 }' sjoint.txt > chr9.teos.txt

awk '$2 == "10" { print $0 }' sjoint.txt > chr10.teos.txt


## Decreasing value 

sort -k3,3nr joint.txt >  dt.txt

awk '$2 == "1" { print $0 }' dt.txt > dchr1.teos.txt 		##chr1

awk '$2 == "2" { print $0 }' dt.txt > dchr2.teos.txt 		##chr2

awk '$2 == "3" { print $0 }' dt.txt > dchr3.teos.txt

awk '$2 == "4" { print $0 }' dt.txt > dchr4.teos.txt

awk '$2 == "5" { print $0 }' dt.txt > dchr5.teos.txt

awk '$2 == "6" { print $0 }' dt.txt > dchr6.teos.txt

awk '$2 == "7" { print $0 }' dt.txt > dchr7.teos.txt

awk '$2 == "8" { print $0 }' dt.txt > dchr8.teos.txt

awk '$2 == "9" { print $0 }' dt.txt > dchr9.teos.txt

awk '$2 == "10" { print $0 }' dt.txt > dchr10.teos.txt


## unknown positions 

awk '$3 == "unknown" { print $0 }' dt.txt > uk.teos.txt



## move all teosinte files to a new folder 

mkdir teosfolder

mv chr1.teos.txt chr2.teos.txt chr3.teos.txt chr4.teos.txt chr5.teos.txt chr6.teos.txt chr7.teos.txt chr8.teos.txt chr9.teos.txt chr10.teos.txt teosfolder/

mv dchr1.teos.txt dchr2.teos.txt dchr3.teos.txt dchr4.teos.txt dchr5.teos.txt dchr6.teos.txt dchr7.teos.txt dchr8.teos.txt dchr9.teos.txt dchr10.teos.txt uk.teos.txt teosfolder/


```









	

 

