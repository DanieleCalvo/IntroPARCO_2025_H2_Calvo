# IntroPARCO_2025_H2_Calvo
This project has the purpose to undertand if it's better to parallelize the project instead of keeping it sequential and what parallelization logic is better, knowing the number of threads available in the hardware.

# What to do
## Installing GlobalProtect
[Click here to install the VPN necessary to use the cluster safely](https://servicedesk.unitn.it/sd/it/kb-article/globalprotect-palo-alto-client-downloads?id=kb_article_view&sysparm_article=KB0011269)

when you open it, put as portal `vpn-mfa.icts.unitn.it`

## Installing MobaXTerm
MAC and Linux provide a built-in SSH client, but I'm on windows so I installed Mobaxterm (PuTTY is a valid alternative). To install mobaxterm, visit [the official web page](https://mobaxterm.mobatek.net/)
## Accessing the cluster
with the terminal (for MAC and Linux) or mobaxterm open, digit `ssh YourName.YourSurname@hpc.unitn.it` and with the password of the istitutional mail you can access the cluster

<img src="Images_readme/ssh.jpg" alt="ssh" width="500">

## Compiling the `.pbs` file
Move the desired directory into your mobaxterm

<img src="Images_readme/move.jpg" alt="move" width="500">

Then, digit into the terminal cd `Name of the directory you want to use`

<img src="Images_readme/cd.jpg" alt="cd" width="500">

use ls to see the files available

<img src="Images_readme/ls.jpg" alt="ls" width="500">

after you saw the name of the pbs file, open it clicking it on the left (where you moved it), so you can change the path to  `/home/YourName.YourSurname/DirectoryName/`
then, use qsub `NameOfThePbsFile.pbs` to run it

<img src="Images_readme/qsub.jpg" alt="qsub" width="500">

with qstat -u `YourName.YourSurname` you can see the status of your pbs

<img src="Images_readme/qstat.jpg" alt="qstat" width="500">

Sometime the pbs refuse to collaborate. It could be that the file is considered in windows mode. To change it, open it with the default text editor and change it
on top left (you'll see many icons, click on the penguin) (remember to save the file before you close it).

P.S. my results is the merging of more pbs outputs, to have more values and more accurate results
