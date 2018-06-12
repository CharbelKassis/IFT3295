from fasta import Fasta
from FonctionChevauchement import CalculeChevauchement
import numpy as np
import math
import argparse


#constantes
lambd = 0.192
ln_k = math.log(0.176,math.e)
ln_2 = math.log(2,math.e)
match = 5
mismatch = -4
moinsInfini = -999999


def bitscore(score):

    return round( (lambd*score - ln_k) / ln_2)

def eValue(m,n,bitscore):

    return m * n * math.pow(2,-bitscore)

def calculIdent(etend_glob):
    
    for i in range(len(etend_glob)):

        count = 0.0
        mot1 = etend_glob[i]["seq1"]
        mot2 = etend_glob[i]["seq2"]
        
        for j in range(len(mot1)):

            if( mot1[j] == mot2[j] ):

                count = count+1

        etend_glob[i]["ident"] = (count/len(mot1)) * 100

        


#imprime le resultat
def myPrint(etend_glob,fasta,mot1):

    etend_glob = sorted(etend_glob, key=lambda x: x["score"], reverse=True)
    
    for hsp in etend_glob:
    
        d1 = hsp["debut_mot1"]
        f1 = hsp["fin_mot1"]
        d2 = hsp["debut_mot2"]
        f2 = hsp["fin_mot2"]
        
        str_d1 = "{num:02d}".format(num=d1)
        str_f1 = "{num:02d}".format(num=f1)
        str_d2 = "{num:02d}".format(num=d2)
        str_f2 = "{num:02d}".format(num=f2)
        
        seqdb = fasta.getSequence(hsp["mot2_index"])

        print(">" + fasta.getInfo(hsp["mot2_index"]) + " score: "+str(hsp["score_SW"]) + " ident: "+str(hsp["ident"]))
        print(hsp["seq1"])
        print(hsp["seq2"])
        print()
        print("# Best HSP:") 
        print("id:"+fasta.getInfo(hsp["mot2_index"]) + ", score: " + str(hsp["score"]) + " bitscore: " + str(hsp["bitscore"]) + ", evalue: " + str(hsp["evalue"]))
        print(str_d1+" "+mot1[d1:f1]+" "+str_f1)
        print(str_d2+" "+seqdb[d2:f2]+" "+str_f2)
        print()
        print("-"*40)

    print("total : "+str(len(etend_glob)))


def SmithWaterman(etend_glob,mot1,fasta):

    for hsp in etend_glob:
        
        mot2 = fasta.getSequence(hsp["mot2_index"])
        score,sequence2,sequence1 = CalculeChevauchement(mot2,mot1)
        hsp["score_SW"] = score
        hsp["seq1"] = sequence1
        hsp["seq2"] = sequence2
    
def addEValue(etend_glob,m,n):

    for i in range(len(etend_glob)):

        if(etend_glob[i]):
            
            bit_score = etend_glob[i]["bitscore"] = bitscore(etend_glob[i]["score"])
            etend_glob[i]["evalue"] = eValue(m,n,bit_score)

#enlever les hsp qui ont un evalue plus grand que ss
def removeNonPertinant(etend_glob):

    i=0
    while(i<len(etend_glob)):

            if(not etend_glob[i] or etend_glob[i]["evalue"] > ss):
            
                del etend_glob[i]
            else:
                i = i+1

        
#choisir le plus grand score parmis tous les kmer d'un mot en particulier
def reduire(etend_mot):

    if(etend_mot):

        return max(etend_mot, key = lambda x: x["score"])

    return etend_mot

def fusion(resultats,mot1,mot2):

    #mettre ensemble tous les etendements (qu'ils proviennent du meme kmer ou pas)
    all_etendement = []
    
    #pour chaque etendement de kmer different
    for i in range(len(resultats)):

        #pour chaque etendement de kmer similaire
        for j in range(len(resultats[i])):
 
            if(resultats[i][j] != []):
 
                all_etendement.append(resultats[i][j])
    i=0
    #debuter la fusion
    while(i<len(all_etendement)):

        j=1
        while(j<len(all_etendement)):

            if(i != j):

                HSP1_db_debut = all_etendement[i]["debut_mot2"]
                HSP1_db_fin = all_etendement[i]["fin_mot2"]
                HSP2_db_debut = all_etendement[j]["debut_mot2"]
                HSP2_db_fin = all_etendement[j]["fin_mot2"]
                HSP1_input_debut = all_etendement[i]["debut_mot1"]
                HSP1_input_fin = all_etendement[i]["fin_mot1"]
                HSP2_input_debut = all_etendement[j]["debut_mot1"]
                HSP2_input_fin = all_etendement[j]["fin_mot1"]

                cond1= HSP1_db_fin >= HSP2_db_debut or HSP1_db_debut <= HSP2_db_fin
                cond2= HSP1_input_fin >= HSP2_input_debut or HSP1_input_debut <= HSP2_input_fin
                cond3= HSP2_db_debut - HSP1_db_debut == HSP2_input_debut - HSP1_input_debut

                #si les 3 conditions sont satisfaites, alors on peut fusioner les deux kmer etendus
                if(cond1 and cond2 and cond3):

                    all_etendement[i]["debut_mot1"] = min(HSP1_input_debut,HSP2_input_debut)
                    all_etendement[i]["fin_mot1"] = max(HSP1_input_fin,HSP2_input_fin)
                    all_etendement[i]["debut_mot2"] = min(HSP1_db_debut,HSP2_db_debut)
                    all_etendement[i]["fin_mot2"] = max(HSP1_db_fin,HSP2_db_fin)
                    del all_etendement[j]
                    
                #s'il n'y a pas de fusion, alors voir le prochain etendement    
                else:
                    j = j+1
                    
            #si i=j, aller verifier le prochain etendement
            else:
                j = j+1     

        i = i+1

    #recalculer score
    for i in range (len(all_etendement)):

        deb_mot1 = all_etendement[i]["debut_mot1"]
        fin_mot1 = all_etendement[i]["fin_mot1"]
        deb_mot2 = all_etendement[i]["debut_mot2"]
        fin_mot2 = all_etendement[i]["fin_mot2"]
        fin = all_etendement[i]["fin_mot1"] - deb_mot1 + 1

        j = 0
        score = 0
        while(j < fin):

            score = score + delta(mot1[deb_mot1+j],mot2[deb_mot2+j])
            j = j+1

        all_etendement[i]["score"] = score

    return all_etendement

    
#fonction servant a entendre un kmer dans mot2
#kmer = le kmer
#kmer_pos = ou commence le kmer dans mot1
#mot1 = le mot duquel on a extrait le kmer
#mot2 = le mot dans lequel on cherche un match avec le kmer
#mot2_index = l'endroit ou le mot se trouve dans la db
#seed = seed
#E = seuil
def etendre(kmer,kmer_pos,mot1,mot2,mot2_index,seed,E):
 
    index = index_of(kmer,mot2,seed) #trouver la ou un kmer match dans mot2 en prenant en compte le seed

    results = [] 
    #le resultat de chaque etendement. Chaque entree est de cette forme : {kmer,debut_mot1,fin_mot1,debut_mot2,fin_mot2,mot2_index,score}
    #kmer = la ou se trouve le kmer dans mot1
    #debut_mot1 = le debut du resultat de l'etendement dans mot1
    #fin_mot1 = la fin du resultat de l'etendement dans mot1
    #debut_mot2 = le debut du resultat de l'etendement dans mot2
    #fin_mot2 = la fin du resultat de l'etendement dans mot2
    #mot2_index = la position du mot2 dans la db
    #score = le score

    if(index):
        
        for i in range (len(index)):

            debut_mot1 = kmer_pos
            fin_mot1 = kmer_pos+len(kmer)-1
            debut_mot2 = index[i]
            fin_mot2 = index[i]+len(kmer)-1
            max_score = len(kmer)*match
            current_score = max_score
            
            while(True):

                score_d = scoreDroite(mot1,mot2,fin_mot1,fin_mot2)
                
                if(not isScoreAcceptable(current_score,score_d,max_score,E)): #On n'a pas depasse le seuil par la droite
                    score_d = moinsInfini
                    
                score_g = scoreGauche(mot1,mot2,debut_mot1,debut_mot2)

                if(not isScoreAcceptable(current_score,score_g,max_score,E)): #On n'a pas depasse le seuil par la gauche
                    score_g = moinsInfini

                score_dg = score_d + score_g
                
                #il n'y a eu aucun etendement, donc l'etendement est termine pour ce kmer
                if(score_d == moinsInfini and score_g == moinsInfini): 
                    
                    results.append({"kmer":kmer_pos, "mot2_index":mot2_index, "debut_mot1":debut_mot1, "fin_mot1":fin_mot1, "debut_mot2":debut_mot2, "fin_mot2":fin_mot2, "score":max_score})
                    break #sortir de la boucle pour aller etendre le prochain kmer, s'il y en a
                
                else:
                   
                    score_array = np.array([score_g,score_d,score_dg]) #mettre les 3 scores dans un tableau numpy)
                    max_value = np.amax(score_array) #aller chercher le maximum
                    index_max = np.argmax(score_array) #ainsi que l'index du maximum
                    current_score = current_score + max_value #MAJ de current score
                    max_score = max(max_score,current_score) #MAJ de max score

                    if(index_max == 0): #si score max vient de gauche, etendre a droite

                        debut_mot1 = debut_mot1 - 1
                        debut_mot2 = debut_mot2 - 1

                    elif(index_max == 1): #sinon si score vient de droite, etendre gauche

                        fin_mot1 = fin_mot1 + 1
                        fin_mot2 = fin_mot2 + 1  

                    else: #sinon etendre des deux cotes

                        fin_mot1 = fin_mot1 + 1
                        fin_mot2 = fin_mot2 + 1
                        debut_mot1 = debut_mot1 - 1
                        debut_mot2 = debut_mot2 - 1

    return results

#verifier si on a depasser le seuil
def isScoreAcceptable(current_score,new_score,max_score,E):

    return current_score + new_score - max_score > -E
            

#calculer le score avec un etendement a droite
def scoreDroite(mot1,mot2,fin_mot1,fin_mot2) :
    
    #verifier si on ne deborde de mot1 ou de mot2 par la droite
    if(fin_mot1+1 != len(mot1) and fin_mot2+1 != len(mot2)):

        return delta(mot1[fin_mot1 + 1],mot2[fin_mot2 + 1])
        
    return moinsInfini
           

#calculer le score avec un etendement a gauche
def scoreGauche(mot1,mot2,debut_mot1,debut_mot2):

    #verifie si on ne deborde pas de mot1 ou de mot2 par la gauche
    if(debut_mot1 != 0 and debut_mot2 != 0):
        
        return delta(mot1[debut_mot1 - 1],mot2[debut_mot2 - 1])
    
    return moinsInfini

#retourne le score de match si v_i = w_j, le score de mismatch sinon
def delta(v_i,w_j):

    if(v_i == w_j):

        return match

    return mismatch


#trouver un match du mot1 dans mot2 en prenant compte du seed, algorithme naive o(n^2)
def index_of(mot1,mot2,seed):

    if(len(seed) != len(mot1)): #il faut que le seed soit de la meme taille que le mot
        return None 

    index = [] #contient les index ou il y a un match

    for i in range (len(mot2)-len(mot1)+1):

        match = True
        
        for j in range (len(mot1)):
            if(seed[j] == "1" and mot2[i+j] != mot1[j]):
                match = False
                break
        if(match):
            index.append(i)

    return index

#genere tous k-mer de la sequence seq
def getKMers(seq,k):

    kmers = []
    
    for i in range(0,len(seq)-k+1):
        kmers.append(seq[i:i+k])

    return kmers


#main
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",help="la sequence en input",required=True)
    parser.add_argument("-db",help="chemin(relatif ou absolu) vers le fichier fasta",required=True)
    parser.add_argument("-E",help="seuil pour le score lors d'un etendement",type=float,default=4)
    parser.add_argument("-ss",help="seuil de pertinance pour le evalue",type=float,default=0.001)
    parser.add_argument("-seed","--seed",help="chaine contenant que des 0 ou 1",default="1"*11)
    args = parser.parse_args()
    
    fasta = Fasta(args.db) #charger la db en memoire
    inp = args.i
    m = fasta.charCount()
    n = len(inp)
    E = args.E
    ss = args.ss
    seed = args.seed

    kmers = getKMers(inp,len(seed))

    etend_mot = [] #tableau contiendra les hsp d'un mot particulier dans la db avec tous les kmer
    etend_glob = [] #tableau contiendra tous les hsp

    #pour chaque mot dans la db
    for i in range(len(fasta)):

        #pour chaque kmer
        for j in range(len(kmers)):

            #le resultat entre l'etendement du kmer j avec le mot i dans la db
            resultat = etendre(kmers[j],j,inp,fasta.getSequence(i),i,seed,E)
            
            if(resultat):
                
                etend_mot.append(resultat)
        
        etend_mot = fusion(etend_mot,inp,fasta.getSequence(i)) #fusioner les kmer
        etend_mot = reduire(etend_mot) #choisir celui qui a le plus grand score
        etend_glob.append(etend_mot) #ajouter le resultat dans le tableau qui contient tous les etendements pour tous les mots
        etend_mot = []

    addEValue(etend_glob,m,n) #ajouter bitscore et evalue a chaque hsp
    removeNonPertinant(etend_glob) #retirer les hsp non pertinants en utilisant le seuil ss
    SmithWaterman(etend_glob,inp,fasta) #ajouter le score de l'alignement local ainsi que les alignements
    calculIdent(etend_glob) #ajouter le pourcentage d'identite de la region de similarite maximale
    
    myPrint(etend_glob,fasta,inp)


        


        
        
    
