
I = [];
Q = [];
mI = [];
mQ = [];
frqBinIndex = 1;
codePhase1 = codePhase;
for i = 1:1:19
A = eval(['CM' num2str(i,'%02d')]);
['CM' num2str(i,'%02d')]
i1 = real(A(frqBinIndex,codePhase1));
q1 = imag(A(frqBinIndex,codePhase1));
% i1 = real(A(fd(i+1),cp(i+1)));
% q1 = imag(A(fd(i+1),cp(i+1)));

I = [I;i1];
Q = [Q;q1];
end
Im = mean(I);
Qm = mean(Q);
% plot(1:21,I)
scatter(I,Q)
hold on
scatter(Im,Qm,'+')
grid on
axis equal
% axis square
xlabel('Inphase component')
ylabel('Quadrature component')
% mI = [mI,Im];
% mQ = [mQ,Qm];