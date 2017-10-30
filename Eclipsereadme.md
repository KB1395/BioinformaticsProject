#FSH - Fast Spaced-seed Hashing

##Description

Patterns with wildcards in specified positions, namely spaced seeds, are increasingly used instead of k-mers in many bioinformatics applications that require indexing, querying and rapid similarity search, as they can provide better sensitivity. Many of these applications require to compute the encoding or the hashing of each position in the input sequences with respect to the given spaced seed, or to multiple spaced seeds. While the hashing of k-mers can be rapidly computed by exploiting the large overlap between consecutive k-mers, spaced seeds hashing is usually computed from scratch for each of the positions in the input sequence, thus resulting in slower processing.  
  
Fast Spaced-seed Hashing (FSH), exploits the similarity of the hash values of spaced seeds computed at adjacent positions in the input sequence, similarly to rolling hash. We also propose a generalized version of the algorithm for the simultaneous computation of multiple spaced seeds hashing.  
  
In our experiments, FSH can compute the hashing values of spaced seeds with a speedup,  with respect to the traditional approach, between 1.6x to 5.3x, depending on the structure of the spaced seed.  
  
Spaced seed hashing is a routine task for several bioinformatics application. FSH allows to perform this task efficiently. This has the potential of major impact in the field, making spaced seed applications not only accurate, but also faster and more efficient. 

#Option Project Eclipse
Link library if you need:  
+ GCC C++ Linker -> miscellaneus -> -fopenmp  
  
CPP compiler options: -std=c++0x -O3 -Wall -c -fmessage-length=0 (-fopenmp optional)  
Dialect: -std=c++0x  

#Download
*NB:This software tests only the computation time and give the results in files*

##Compilation:
Open terminal and go to FSH/Release/ and then:  
make all  

###Option compilation
Y = Number of core  
-jY

#Algorithm Option
In the following paragraph is described the input file's structure and the parameters available in the FSH algorithm.  

##File accepted and structure
File accepted have the following structure:  
Structure file .fna example:  
> \>IDENTIFICATION  
> ATAATTGGCAAGTGTTTTAGTCTTAGAGAGATTCTCTAAGTCTAACTTGAACTCAATTTGGAATCATTTCCCAATTTTTA

Structure .fastq exemple:  
> @IDENTIFICATION #0/1  
> CCCATGCCTTTAGCCAAATTCACGGTTTGATCACCCCTAAAACCAGCCAATATACCGAAGTGGAAGCCAGCATAAATGGCCTCAATATTACCGAAATGGAT  
> +  
> HBIIIIIIIHHDIHIGIIGGIHIIGIDIIIIBIHI@IIH@HIIHIIF5IIHEII>BDAHIBIEDBEIDG@HAEH*I@AEI=#CE?G17EEDHDEB@@?#8B  
  
In this NEW VERSION, paired-end reads are passed to the algorithm in separeted file in which the reads are paired 
in the same order in which are writen. So we raccomend to control the paired-end read if they are paired in the
correct manner.

##Parameter
**-si** File path single-end reads  
**-pi** File paths paired-end reads   
**-dirO** Path directory output files. Default: output/  
**-q** Enter a spaced seeds path as -q <AbsPathFile>. Every file's line must contain a spaced qmer. Spaced qmer used for overlapping reads. Ex. 1\*\*\*1\*111. 1 is the simbol considered, any others are not valid. Default spaced seeds are:  
_1111011101110010111001011011111_ -> CLARK-S paper  
_1111101011100101101110011011111_ -> CLARK-S paper  
_1111101001110101101100111011111_ -> CLARK-S paper  
_1111010111010011001110111110111_ -> rasbhari minimizing overlap complexity  
_1110111011101111010010110011111_ -> rasbhari minimizing overlap complexity  
_1111101001011100111110101101111_ -> rasbhari minimizing overlap complexity  
_1111011110011010111110101011011_ -> rasbhari maximizing sensitivity  
_1110101011101100110100111111111_ -> rasbhari maximizing sensitivity  
_1111110101101011100111011001111_ -> rasbhari maximizing sensitivity  

#Run
Calls algorithm where is compiled:  
./FSH -si ../TestInputFile/long_example_1.fna -q 11101110110110111101
./FSH -pi ../TestInputFile/short_example_1.fna.1 ../TestInputFile/short_example_1.fna.2 -q 11101110110110111101  
./FSH -pi ../TestInputFile/short_example_2.fna.1 ../TestInputFile/short_example_2.fna.2 -q 11101110110110111101 .dirO /home/user/desktop/  