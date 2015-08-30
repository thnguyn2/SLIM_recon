%4-digit file Read
function A = RawReadREPZTimelapse(fname,imRep, num,znum, channel)

if imRep<10; numlength=4;
elseif (imRep>=10 & imRep<100); numlength=3;
elseif (imRep>=100 & imRep<1000); numlength=2;
elseif (imRep>=1000 & imRep<10000); numlength=1;
else numlength=0;
end
for m=1:numlength;
     mm=num2str(0);
     fname=strcat(fname,mm);
end
mm=num2str(imRep);
fname=strcat(fname,mm,'t0000');

if num<10; numlength=3;
elseif (num>=10 & num<100); numlength=2;
elseif (num>=100 & num<1000); numlength=1;
else (num>=1000 & num<10000); numlength=0;
end

for m=1:numlength;
     mm=num2str(0);
     fname=strcat(fname,mm);
end
mm=num2str(num);
fname=strcat(fname,mm,'z');

if znum<10; numlength=2;
elseif (znum>=10 & znum<100); numlength=1;
elseif (znum>=100 & znum<10000); numlength=0;
else (znum>=1000 & znum<10000); numlength=0;
end

for m=1:numlength;
     mm=num2str(0);
     fname=strcat(fname,mm);
end

     

     mm=num2str(znum);
     fname=strcat(fname,mm,channel,'.tif');
     A=im2double(imread(fname));
return;
