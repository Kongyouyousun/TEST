
SingleOFDM=var(pdschGrid(:,1:12));
alpha=2048/(512+2048);
energyfft=0;
energy=0;
for i=1:12
    energyfft=energyfft+alpha*SingleOFDM(:,i);
    energy=energy+SingleOFDM(:,i);
end  
energyfft
energy

nablda=energyfft/energy














