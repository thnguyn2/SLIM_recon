function a=TrRead(filename)
%fid=fopen('abc.txt');
fid=fopen(filename);
m=1;
n=1;
while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end
  % disp(tline)
    if isempty(tline)
        tline=fgetl(fid);
        m=m+1
        n=1;
    elseif tline(1)~='%'
    a{m}(n,:)=str2num(tline);
    n=n+1;
    end
    
end
fclose(fid)