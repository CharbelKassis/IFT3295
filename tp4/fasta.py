#Classe servant a acceder aux sequences dans un fichier Fasta

indel = "-"

class Fasta:

    def __init__(self,filename):

        self.filename = filename
        self.__info = [] #contiendra les informations sur la sequence i
        self.__sequences = [] #contiendra la sequence i
        self.__map = {} #map {nom_sequence : sequence}
        self._getFasta()
      
    def _getFasta(self):
  
        fastaFile = open(self.filename,"r")

        sequence = ""

        for line in fastaFile:

            if(line[0] == ">"): #information sur la sequence courante

                if(sequence): #verifier qu'on a pas commencer l'algorithme pour ajouter la sequence courrante dans la liste des sequences

                    self.__sequences.append(sequence)
                    sequence = ""
                    
                line = line[1:-1] #enlever > et le \n
                self.__info.append(line)

            else:
                
                line = line[:-1] #enlever le caracter \n 
                sequence = sequence + line
             
        self.__sequences.append(sequence) #ajouter la derniere sequence
        fastaFile.close()

        for i in range(len(self.__sequences)):

            self.__map[self.__info[i]] = self.__sequences[i]
             
    def info(self,i):
 
        return self.__info[i]

    #fonction servant a enlever les positions des sequences ou il y a un gap
    def filtrer(self):

        indel_pos = set()
        
        seq_len = self.seq_len()

        #pour chaque sequence
        for seq in self:
        
            #pour chaque position de sequence
            for i in range(seq_len):

                #s'il y a un indel
                if(self[seq][i] == indel):

                    #mettre la position en memoire
                    indel_pos.add(i)

        #pour chaque sequence
        for seq in self:

            newSequence = ""
            
            #pour chaque position de sequence
            for i in range(seq_len):

                #si la position n'a pas de gap
                if(not i in indel_pos):

                    #l'ajouter a newSequence
                    newSequence += self[seq][i]

            self[seq] = newSequence

    #retourne la taille des sequences
    def seq_len(self):

        return len(self.__sequences[0])

    def __len__(self):

        return len(self.__sequences)

    def __iter__(self):

        return iter(self.__info)

    def __getitem__(self,key):
 
        return self.__map[key]

    def __setitem__(self,key,value):

        self.__map[key] = value

    

    
