%4-digit file Read
function A = RawRead(fname, num)

if num<10; numlength=2;
elseif (num>=10 & num<100); numlength=1;
elseif (num>=100 & num<10000); numlength=0;
else (num>=1000 & num<10000); numlength=0;
end

for m=1:numlength;
     mm=num2str(0);
     fname=strcat(fname,mm);
end
     mm=num2str(num);
     fname=strcat(fname,mm,'.tif');
     A=im2double(imread(fname));
return;

