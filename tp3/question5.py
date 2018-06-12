from fasta import Fasta
import subprocess
import os

unknown = Fasta("unknown.fasta")
db = "tRNA.fasta"
cur_dir = os.getcwd()

if not os.path.exists(os.path.join(cur_dir,"question5")):
    os.makedirs(os.path.join(cur_dir,"question5"))

print("Question 5")
print()
    
for i in range(len(unknown)):
    
    path = os.path.join("question5","output"+str(i)+".txt")
    inp = unknown.getSequence(i)
    process = subprocess.Popen("python plast.py -i "+inp+" -seed 111010010100110111"+" -ss 0.01"+" -db "+db+" >> "+path , shell=True)
    process.wait()
    print("calcul pour "+inp+" termine. Le resultat est dans le fichier "+path)