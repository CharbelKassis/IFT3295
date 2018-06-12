#classe servant a stocker les donnes Fasta en memoire afin de les utiliser dans un autre programme

class Fasta:

    def __init__(self,filename):
        self.filename = filename
        self.__info = [] #contiendra les informations sur la sequence i
        self.__sequences = [] #contiendra la sequence i
        self.__nbOfChar = 0 #nombre de caracteres total de toute la base de donne
        self._getFasta()
        
        
    def _getFasta(self):
 
        seqFile = open(self.filename,"r")

        for line in seqFile:
        
            if(line[0] == ">"): #cette ligne indique l'information sur la sequence
                
                line = line[1:-1] #enlever > et le \n
                self.__info.append(line)
                
            else:

                line = line[:-1] #enlever \n
                self.__sequences.append(line)
                self.__nbOfChar = self.__nbOfChar + len(line) #maj le nombre de caracteres
                
        seqFile.close()
        
    def getSequence(self,i):

        return self.__sequences[i]

    def getInfo(self,i):

        return self.__info[i]

    def __len__(self):

        return len(self.__info)

    def charCount(self):

        return self.__nbOfChar

     
          
               
