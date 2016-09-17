%This is the main script for the generation of a new duplex.
%It requires as INPUT, Z,T01, T10, T11, D01b, D10a, D11a, D11b, matrices
%that are the results of "code_entropy_duplex.m"
% P00, P01, P10, P11 are matrices that encode the multilink probabilities,
% for each link
%they are transformed in vectors e joined together in a single matrix p
% This script calls the function "duplexsingleinstance.m"
%OUTPUT the variable duplex stores in the first column Layer 1 and in the
%second column Layer 2, the adjacency matrix for both layers can be
%obtained using the command squareform

P01=(T01./Z)./(D01b+(D01b==0));
P01=P01-diag(diag(P01));
P01=squareform(P01)';

P10=(T10./Z)./(D10a+(D10a==0));
P10=P10-diag(diag(P10));
P10=squareform(P10)';

P11=(T11./Z)./(D11a+(D11a==0))./(D11b+(D11b==0));
P11=P11-diag(diag(P11));
P11=squareform(P11)';

P00=1-P11-P10-P01;

p=[P00,P10, P01, P11];

D01b=D01b-diag(diag(D01b));
d01b=squareform(D01b);

D10a=D10a-diag(diag(D10a));
d10a=squareform(D10a);

D11a=D11a-diag(diag(D11a));
d11a=squareform(D11a);

D11b=D11b-diag(diag(D11b));
d11b=squareform(D11b);


duplex=duplexsingleinstance(p, d10a, d01b, d11a, d11b);