from EdgeList import calculTabReads

#Code ecrit en python36
#Par Charbel Kassis et Boumediene Boukharouba
#TP1
#Question 2

#Ce programme produira le fichier edgelist.txt a partir de reads.fq
input_file = open('reads.fq','r')
output_file = open('edgelist.txt','w')
calculTabReads(input_file,output_file)


