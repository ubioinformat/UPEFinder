# UPEfinder

The Human Proteome Project (HPP) is leading the international effort to characterize the human proteome. Although the main goal of this
project was first focused on the detection of the missing proteins, a new challenge arose from the need to assign biological functions to 
the uncharacterized human proteins and describe their implication in human diseases. Not only the proteins with experimental evidence 
(uPE1s) were the object of study of this challenge, neXt-CP50, but also the uncharacterized missing proteins (uMPs). In this work we 
developed a new bioinformatic approach to infer biological annotations for the uPE1 proteins based on a "guilt-by-association" analysis.
We used the correlation of these proteins with the well characterized PE1 proteins to construct a network. In this way, we applied the 
well-known Page Rank algorithm to this network identifying the most relevant nodes, which were the biological annotations of the uPE1 
proteins. The RNA-Seq datasets analyzed and all the generated information were stored in a database. In addition, we implemented the web
application UPEFinder (https://upefinder.proteored.org) in order to facilitate the access to this new resource. This information is 
specially relevant for the researchers of the HPP project interested in the generation and validation of new hypotheses about the 
functions of the uPE1 proteins. Both the database and the web application are publicly available (GITHUB).
