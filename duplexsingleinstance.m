% See the script "main_single_instance.m" for explanations
function [duplex]=duplexsingleinstance(p, d10a, d01b, d11a, d11b)
l=size(p,1);
duplex=zeros(l,2);


for ind=1:l
    r =find(mnrnd(1,p(ind,:)));
    switch r
        case 1
            duplex(ind,:)=[0 0];
        case 2
            duplex(ind,:)=[1+(geornd(d10a(ind)/(1+d10a(ind)))) 0];
        case 3
            duplex(ind,:)=[0 1+(geornd(d01b(ind)/(1+d01b(ind))))];
        otherwise
            duplex(ind,:)=1+geornd([d11a(ind)/(1+d11a(ind)), d11b(ind)/(1+d11b(ind))]);
    end
    
    
    
end