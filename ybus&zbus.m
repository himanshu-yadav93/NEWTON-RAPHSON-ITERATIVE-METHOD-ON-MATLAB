% program to from admittance bus formation with transformer tap setting
linedata=linedata()
fb=linedata(:,1);
tb=linedata(:,2);
r=linedata(:,3);
x=linedata(:,4);
b=linedata(:,5);
a=linedata(:,6);
z=r+1i*x;
y=1./z;
b=1i*b;
nbus=max(max(fb),max(tb));
nbranch=length(fb);
ybus=zeros(nbus,nbus);
  
for k=1:nbranch 
    ybus(fb(k),tb(k))=ybus(fb(k),tb(k))-y(k)/a(k);
    ybus(tb(k),fb(k))=ybus(fb(k),tb(k));
end

for m=1:nbus
    for n=1:nbranch
        if fb(n)==m
            ybus(m,m)=ybus(m,m)+y(n)/(a(n)^2)+b(n);
        elseif tb(n)==m 
            ybus(m,m)=ybus(m,m)+y(n)+b(n);
        end
    end
end
 zbus=inv(ybus);

 disp('ybus matrix:');
 disp(ybus);
 disp('zbus matrix:');
 disp(zbus);