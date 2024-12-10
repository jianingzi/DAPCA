% A robust PCA algorithm for drift fault detection

function [FAR,F1,fai,kesi,result,AUC]=RPCA(Y,Xtest1,label,omega,k,i)

alpha=k(i);
% training
[X_row,X_col] = size(Y);
sigmaXtrain = cov(Y);
[T,lamda] = eig(sigmaXtrain);
D = flipud(diag(lamda));
num_pc = 1;
while sum(D(1:num_pc))/sum(D) <0.1
    num_pc = num_pc +1;
end
P = T(:,X_col-num_pc+1:X_col);
T2UCL1=num_pc*(X_row-1)*(X_row+1)*finv(alpha,num_pc,X_row - num_pc)/(X_row*(X_row - num_pc));
for i = 1:3
    %        theta(i) = sum((D(num_pc+1:X_col)).^i);
    theta(i) = D(i)/sum(D);
end
h0 = 1 - 2*theta(1)*theta(3)/(3*theta(2)^2);
ca = norminv(alpha,0,1);
QUCL = theta(1)*(h0*ca*sqrt(2*theta(2))/theta(1) + 1 + theta(2)*h0*(h0 - 1)/theta(1)^2)^(1/h0);

% testing
n = size(Xtest1,1);
X_sma = movmean(Xtest1, omega);
[r,y] = size(P*P');
I = eye(r,y);
T2 = zeros(n,1);
Q = zeros(n,1);
fai= zeros(n,1);

for i = 1:n
    T2(i)=X_sma(i,:)*P*pinv(lamda(X_col-num_pc+1:X_col,X_col-num_pc+1:X_col))*P'*X_sma(i,:)';
    Q(i) = norm(X_sma(i,:)*(I - P*P'))^2;
    fai(i)=(Q(i)/QUCL^2)+(T2(i)^2/T2UCL1^2);
end

S=lamda(X_col-num_pc+1:X_col,X_col-num_pc+1:X_col);
FAI=P*pinv(S)*P'/T2UCL1+(eye(X_col)-P*P')/QUCL;
S=cov(Y);
g=trace((S*FAI)^2)/trace(S*FAI);
h=(trace(S*FAI))^2/trace((S*FAI)^2);
kesi =g*chi2inv(alpha,h);

result=zeros(n,1);
[idx]=find(fai>kesi);
result(idx)=1;

% faulty sensor location
for i= 1:n
    C1(i,:)=(X_sma(i,:)*(eye(X_col)-P*P')).^2;
end

C=sum(C1,1);
cr=C/sum(C);
[~,fs]=max(cr);
disp(['the number of the faulty sensor is #',num2str(fs)]);

% output result
TP = sum(result == 1 & label == 1);
FP = sum(result == 1 & label == 0);
TN = sum(result == 0 & label == 0);
FN = sum(result == 0 & label == 1);
FAR = ((FP)/(TN+FP))*100;
p=TP/(TP+FP);
r=TP/(TP+FN);
F1 = 2*(r*p)/(r+p)*100;
[~,~,~,AUC]=perfcurve(label,fai,1);
