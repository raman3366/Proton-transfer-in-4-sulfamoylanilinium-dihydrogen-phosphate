%Chk=CA2_small_opt.chk
%Mem=3000MB
%NProcShared=2
# B3LYP/6-31G* Opt=CalcAll Geom=AllCheckpoint Guess=Read
Integral=(Grid=UltraFine) IOp(7/33=1)
 
 
