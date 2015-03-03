## knowledge-survival

This repo serves as a version control and for my dissertation research, and a great way to inspect my work in progress. The README will also be a makeshift blog where I will post developments to the project. Please feel free to send me an email at bambrose@ucla.edu if you have questions.

####~ 3/2/2015 Learning how to scale a record linking problem ~

In my research I face the issue of finding citation code variations in the Web of Knowledge database. For a haystack of millions of citations there may be several thousand needles, or sets of citations that are probably variable addresses to a single reference. I refer ungenerously to these as coding errors, but from a data cleaning standpoint this is fair, even if these are not always due to transcription errors on the part of WOK coders. To start off my search for a good solution, I am consulting William W. Cohen's "mini-course on record linkage and matching", which with some scrolling can be found [here][miniRL].

A handful from the haystack may look like this:

>...
>FACT ACT INQ BOAR, 1893, 1 PROGR REP FACT ACT, P23    
>FACT ACT INQ BOAR, 1894, 2 PROGR REP FACT ACT, P5     
>FACT INSP COMM PE, 1902, 13 FACT INSP COMM PE, P387   
>FACT INV COMM, 1914, 3 REP FACT INV COMM, P304        
>FAIRP FINN NAT CH, 1924, P FAIRP FINN NAT CHU         
>FAIRP SUOM SYN CH, 1916, P FAIRP SUOM SYN CHU         
>FAIRP SUOM SYN CH, 1926, P FAIRP SUOM SYN CHU         
>FAIRP SUOM SYN SY, 1919, P FAIRP SUOM SYN SYN         
>FAIRP SUOM SYN SY, 1922, P FAIRP SUOM SYN SYN         
>FAM WELF ASS AM, 1929, PREL REP COMM UNPUB, P39       
>FAM WELF ASS AM, COMM FUT PROGR, P33                  
>FAM WELF ASS AM, DIV WORK PUBL PRIV A                 
>FAM WELF ASS, 1933, UN REL EXP                        
>FARM BOARD, STOKD W, P18                              
>FARM CRED ADM, 1934, 2 FARM CRED ADM, P6              
>FARM CRED ADM, 1935, MONTHL REP LOANS DIS             
>...

And a needle would look like this:

>
>DIMAGGIO P., 1983, AM SOCIOL REV, V47, P147
>DIMAGGIO P.J., 1983, AM SOCIOL REV, V48, P47
>DIMAGGIO PJ, 1983, AM SOCIOL REV, V48, P147
>

My first attempt to solve the problem of finding these variation sets was inefficient. What was convenient to program in R using the [stringdist package][stringdist] was very computationally inefficient. Most obviously, I didn't take advantage of the fact that pairs are unordered for some string measures like Jaro-Winkler, and my implementation calculated the distance twice, once for each ordered pair. That embarassing misstep is easily avoided, but there are even more sophisticated approaches that eliminate redundancies even between pairs. These approaches rely on sorting the list and carrying results forward when the calculation is identical in a subsequent pair, as when the initial substrings of a series of strings are identical. A data structure called a trie is the solution, as explained [here][trie1].

I hope to have something in the works soon, and to put it up on [Savio][savio] shortly!

####~ 2/13/2015 [D-Lab Social Computing Working Group][SCWG] ~

This talk provided an introduction to my coding workflow and illustrated some of the data management and research design decisions that face anyone who studies sociocultural networks. Take a look at the [prezi][2-13-15SCWG] for a visual guide to my suite of programming functions. I hope to add to this prezi as a way of explaining the entire workflow of the project.

[2-13-15SCWG]:https://prezi.com/hf-1-nca8kky/wok/
[SCWG]:http://dlab.berkeley.edu/working-groups/social-computing
[miniRL]:http://www.cs.cmu.edu/~wcohen/
[stringdist]:https://github.com/markvanderloo/stringdist
[trie1]:http://stevehanov.ca/blog/index.php?id=114
[savio]:http://research-it.berkeley.edu/services/high-performance-computing