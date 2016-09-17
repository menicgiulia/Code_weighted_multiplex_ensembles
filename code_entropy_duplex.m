%Matlab function for the comuputation of the entropy value and the related
%Lagrangian multipliers for a canonical duplex ensemble with given expected
%multidegree and multistrength
%INPUT
%multidegree k01, k10, k11, each of them is a column vectors (number of
%nodes x 1)
%multistrength s01b, s10a, s11a, s11b, each of them is a column vectors (number of
%nodes x 1)
%01, 10, 11 are associated with the possible multilinks while a,b are
%associated to Layer 1 and Layer 2
%OUTPUT
%S Entropy value
%t01, t10, t11 Lagrangians multipliers releted to the multidegree sequence
%these are functions of the real Lagrangian multipliers (y=e^(-x) with x the associated true Lagrangian multiplier)
%z01b_new, z10a_new, z11a_new, z11b_new Lagrangians multipliers releted to the multistrength sequence
%these are functions of the real Lagrangian multipliers (y=e^(-x) with x the associated true Lagrangian multiplier)
%k01cal, k10cal, k11cal, s01bcal, s10acal, s11acal, s11bcal fixed average
%values of the constraints, computed with the Lagrangian multipliers
%Z,T01, T10, T11, D01b, D10a, D11a, D11b matrices releted to the partition
%function and the Lagrangian multipliers, useful to speed up the generation
%of new duplexes, belonging to the ensemble 


function [S,t01, t10, t11, z01b_new, z10a_new, z11a_new, z11b_new, k01cal, k10cal, k11cal, s01bcal, s10acal, s11acal, s11bcal,Z,T01, T10, T11, D01b, D10a, D11a, D11b]=code_entropy_duplex(k01, k10, k11, s01b, s10a, s11a, s11b)

precision=10^(-4);
loops=100000;
n=length(k01);



t01=rand(n,1);
t10=rand(n,1);
t11=rand(n,1);

z01b=rand(n,1);
z10a=rand(n,1);
z11a=rand(n,1);
z11b=rand(n,1);


oldt01=zeros(n,1);
oldt10=zeros(n,1);
oldt11=zeros(n,1);

oldz01b=zeros(n,1);
oldz10a=zeros(n,1);
oldz11a=zeros(n,1);
oldz11b=zeros(n,1);

for kk=1:loops
    display(kk)
    
    
    T01=t01*t01';
    T10=t10*t10';
    T11=t11*t11';
    
    D01b=(z01b*ones(1,n))+ (z01b*ones(1,n))'+ z01b*z01b';
    D10a=(z10a*ones(1,n))+ (z10a*ones(1,n))'+ z10a*z10a';
    D11a=(z11a*ones(1,n))+ (z11a*ones(1,n))'+ z11a*z11a';
    D11b=(z11b*ones(1,n))+ (z11b*ones(1,n))'+ z11b*z11b';
    
    %partition function
    Z=1+ T01./(D01b+(D01b==0)) +   T10./(D10a+(D10a==0)) + T11./(D11a.*D11b+((D11a.*D11b)==0));
    
  
    %p01
    num=ones(n,1)*t01';
    den=Z.*D01b;
    summat=num./(den+(den==0));
    summat=summat-diag(diag(summat));
    summat=sum(summat,2);
    t01=k01./(summat+(summat==0));
    
    T01=t01*t01';
    Z=1+ T01./(D01b+(D01b==0)) +   T10./(D10a+(D10a==0)) + T11./(D11a.*D11b+((D11a.*D11b)==0));

    
    %p10
    num=ones(n,1)*t10';
    den=Z.*D10a;
    summat=num./(den+(den==0));
    summat=summat-diag(diag(summat));
    summat=sum(summat,2);
    t10=k10./(summat+(summat==0)); 
    
    T10=t10*t10';
    Z=1+ T01./(D01b+(D01b==0)) +   T10./(D10a+(D10a==0)) + T11./(D11a.*D11b+((D11a.*D11b)==0));

    %p11
    num=ones(n,1)*t11';
    den=Z.*D11a.*D11b;
    summat=num./(den+(den==0));
    summat=summat-diag(diag(summat));
    summat=sum(summat,2);
    t11=k11./(summat+(summat==0));  
    
    T11=t11*t11';
    Z=1+ T01./(D01b+(D01b==0)) +   T10./(D10a+(D10a==0)) + T11./(D11a.*D11b+((D11a.*D11b)==0));
    
    
    
    %z01b
    num=T01.*(z01b*ones(1,n)).*(1 + (D01b+(D01b==0)).^(-1));
    den=Z.*D01b;
    summat=num./(den+(den==0));
    summat=summat-diag(diag(summat));
    summat=sum(summat,2);
    z01b=summat./(s01b+(s01b==0));
    z01b=max(z01b,10^(-15));
    
    D01b=(z01b*ones(1,n))+ (z01b*ones(1,n))'+ z01b*z01b';
    Z=1+ T01./(D01b+(D01b==0)) +   T10./(D10a+(D10a==0)) + T11./(D11a.*D11b+((D11a.*D11b)==0));

    
    
    %z10a
    num=T10.*(z10a*ones(1,n)).*(1 + (D10a+(D10a==0)).^(-1));
    den=Z.*D10a;
    summat=num./(den+(den==0));
    summat=summat-diag(diag(summat));
    summat=sum(summat,2);
    z10a=summat./(s10a+(s10a==0)); 
    z10a=max(z10a,10^(-15));
    
    D10a=(z10a*ones(1,n))+ (z10a*ones(1,n))'+ z10a*z10a';
    Z=1+ T01./(D01b+(D01b==0)) +   T10./(D10a+(D10a==0)) + T11./(D11a.*D11b+((D11a.*D11b)==0));

    
    
    
    %z11a
    num=T11.*(z11a*ones(1,n)).*(1 + (D11a+(D11a==0)).^(-1));
    den=Z.*D11a.*D11b;
    summat=num./(den+(den==0));
    summat=summat-diag(diag(summat));
    summat=sum(summat,2);
    z11a=summat./(s11a+(s11a==0));    
    z11a=max(z11a,10^(-15));  
    
    D11a=(z11a*ones(1,n))+ (z11a*ones(1,n))'+ z11a*z11a';
    Z=1+ T01./(D01b+(D01b==0)) +   T10./(D10a+(D10a==0)) + T11./(D11a.*D11b+((D11a.*D11b)==0));

    

   %z11b
    num=T11.*(z11b*ones(1,n)).*(1 + (D11b+(D11b==0)).^(-1));
    den=Z.*D11a.*D11b;
    summat=num./(den+(den==0));
    summat=summat-diag(diag(summat));
    summat=sum(summat,2);
    z11b=summat./(s11b+(s11b==0));    
    z11b=max(z11b,10^(-15));

    prec2=max(abs((t10>0).*(1-t10./(oldt10+(oldt10==0)))));
    prec3= max(abs((t01>0).*(1-t01./(oldt01+(oldt01==0)))));
    prec4=max(abs((t11>0).*(1-t11./(oldt11+(oldt11==0)))));
    prec5=max(abs((z10a>0).*(1-z10a./(oldz10a+(oldz10a==0)))));
    prec6=max(abs((z01b>0).*(1-z01b./(oldz01b+(oldz01b==0)))));
    prec7=max(abs((z11a>0).*(1-z11a./(oldz11a+(oldz11a==0)))));
    prec8=max(abs((z11b>0).*(1-z11b./(oldz11b+(oldz11b==0)))));


    control=((prec2 < precision) &&(prec3<precision) && (prec4<precision)  && (prec5<precision) && (prec6<precision) && (prec7<precision) )&&(prec8<precision);
    if (control)
           break
    end
   
   
    oldt01=t01;
    oldt10=t10;
    oldt11=t11;

    oldz01b=z01b;
    oldz10a=z10a;
    oldz11a=z11a;
    oldz11b=z11b;
    
    
    
end

display(kk)


    T01=t01*t01';
    T10=t10*t10';
    T11=t11*t11';
    
    D01b=(z01b*ones(1,n))+ (z01b*ones(1,n))'+ z01b*z01b';
    D10a=(z10a*ones(1,n))+ (z10a*ones(1,n))'+ z10a*z10a';
    D11a=(z11a*ones(1,n))+ (z11a*ones(1,n))'+ z11a*z11a';
    D11b=(z11b*ones(1,n))+ (z11b*ones(1,n))'+ z11b*z11b';
    
    %partition function
    Z=1+ T01./(D01b+(D01b==0)) +   T10./(D10a+(D10a==0)) + T11./(D11a.*D11b+((D11a.*D11b)==0));

    Hdegree=0;
    
    den=Z.*D01b;
    summat=T01./(den+(den==0));
    summat=summat-diag(diag(summat));
    k01cal=sum(summat,2);
    Hdegree= Hdegree - summat.*log(T01+(T01==0));
    
    
    
    den=Z.*D10a;
    summat=T10./(den+(den==0));
    summat=summat-diag(diag(summat));
    k10cal=sum(summat,2);
    Hdegree= Hdegree - summat.*log(T10+(T10==0));
    
    
    
    den=Z.*D11a.*D11b;
    summat=T11./(den+(den==0));
    summat=summat-diag(diag(summat));
    k11cal=sum(summat,2);
    Hdegree= Hdegree - summat.*log(T11+(T11==0));
  
    
    Hw=0;
    z01b_new=(1+z01b).^(-1);
    num=T01.*(1 + (D01b+(D01b==0)).^(-1));
    den=Z.*D01b;
    summat=num./(den+(den==0));
    summat=summat-diag(diag(summat));
    s01bcal=sum(summat,2);
    Hw=Hw-log(z01b_new*z01b_new'+((z01b_new*z01b_new')==0)).*summat;


    z10a_new=(1+z10a).^(-1);
    num=T10.*(1 + (D10a+(D10a==0)).^(-1));
    den=Z.*D10a;
    summat=num./(den+(den==0));
    summat=summat-diag(diag(summat));
    s10acal=sum(summat,2);
    Hw=Hw-log(z10a_new*z10a_new'+((z10a_new*z10a_new')==0)).*summat;
 
    
    z11a_new=(1+z11a).^(-1);
    num=T11.*(1 + (D11a+(D11a==0)).^(-1));
    den=Z.*D11a.*D11b;
    summat=num./(den+(den==0));
    summat=summat-diag(diag(summat));
    s11acal=sum(summat,2);
    Hw=Hw-log(z11a_new*z11a_new'+((z11a_new*z11a_new')==0)).*summat;

 
    z11b_new=(1+z11b).^(-1);    
    num=T11.*(1 + (D11b+(D11b==0)).^(-1));
    den=Z.*D11a.*D11b;
    summat=num./(den+(den==0));
    summat=summat-diag(diag(summat));
    s11bcal=sum(summat,2);
    Hw=Hw-log(z11b_new*z11b_new'+((z11b_new*z11b_new')==0)).*summat;

    Sij=Hdegree+Hw+log(Z+(Z==0));
    S=sum(sum(triu(Sij,1)));
    display(S)
    

