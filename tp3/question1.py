from fasta import Fasta
import subprocess
import os

unknown = Fasta("unknown.fasta")
db = "tRNA.fasta"
cur_dir = os.getcwd()

if not os.path.exists(os.path.join(cur_dir,"question1")):
    os.makedirs(os.path.join(cur_dir,"question1"))

print("Question 1")
print()
    
for i in range(len(unknown)):
    
    path = os.path.join("question1","output"+str(i)+".txt")
    inp = unknown.getSequence(i)
    process = subprocess.Popen("python plast.py -i "+inp+" -db "+db+" >> "+path , shell=True)
    process.wait()
    print("calcul pour "+inp+" termine. Le resultat est dans le fichier "+path)
    