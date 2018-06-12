from methodes_utiles import complementReverse
from methodes_utiles import reductionTransitive
from methodes_utiles import writeFASTQ
from methodes_utiles import writeEdge
from edgelist import calculTabReads

#Code ecrit en python36
#Par Charbel Kassis et Boumediene Boukharouba
#TP1
#question 3

#Ce programme produira ces fichiers-la: 
#   le fichier reverse_reads.fq: une version modifiee de reads.fq qui remplace les reads reverses par leur complement-reverse
#	le fichier edgelist_reverse.txt: edgelist genere a partir de reverse_reads.fq
#   le fichier edgelist_reverse_reduite.txt: edgelist des reads sur lesquels on a appliques la reduction transitive 

#tableau qui contient l'index des reads reverse, obtenu dans la question 2
index_reverse = [6,7,9,10,12,15,18]

input_file = open('reads.fq','r')
output_file = open('reverse_reads.fq','w')
input_array = input_file.readlines()

#remplacer les reads reverse par leur complement reverse, sinon reecrire
#ce qui se trouvait deja dans le fichier reads.fq dans un nouveau fichier reverse_reads.fq

print("\nGenerating new FASTQ file reverse_reads.fq")
for i in range(0,len(input_array),4):
    
    sequence = input_array[i+1]

    if (i/4 + 1) in index_reverse:

        sequence = complementReverse(sequence[:-1])+"\n"

    writeFASTQ(output_file,i/4+1,sequence)
print("Finished generating FASTQ file reverse_reads.fq")


input_file.close()
output_file.close()
input_file = open('reverse_reads.fq','r')
output_file = open('edgelist_reverse.txt','w')


scoreMatrix,seqMatrix =  calculTabReads(input_file,output_file)

input_file.close()
output_file.close()

print("\nApplying transitive reduction algorithm on score matrix")
#a partir de la matrice des scores, appliquer l'algorithme de reduction transitive pour obtenir un graphe plus simple
output_file = open('edgelist_reverse_reduite.txt','w')
reductionTransitive(scoreMatrix)

length = len(scoreMatrix)

#ecrire le fichier EdgeList de la matrice des scores sur laquelle la reduction transitive est appliquee
first = True

for i in range(length):
    for j in range(length):
        if(scoreMatrix[i][j] != 0):
            writeEdge(i,j,scoreMatrix,None,output_file)


	



