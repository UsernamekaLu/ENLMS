clear all, clc, close all

N=65;
L=33;
iter=500; %traje dosta vremena
M=10000;
signal=randn(M,1);
signal=signal+sqrt(0.1481)*randn(M,1);
h1=randi([1, 100], 1, 65)';
d=conv(signal,h1);
d=d(1:M);
d=d+sqrt(1e-3)*randn(M,1);

w=zeros(N,1); 
X=zeros(N,L);
e=zeros(M,1);
ksi=zeros(N,1);
z=zeros(N,1);

h=[h1;zeros(N-length(h1),1)];
w_mk_enlms=zeros(M,1);
w_mk_rls=zeros(M,1);
w_mk_nlms=zeros(M,1);
w_mk_lms=zeros(M,1);

for ensemble = 1:iter
    e=zeros(M,1);
    signal=randn(M,1);
    signal=signal+sqrt(0.1481)*randn(M,1);
    d=conv(signal,h1);
    d=d(1:M);
    d=d+sqrt(1e-3)*randn(M,1);
       
    w=zeros(N,1);
    X=zeros(N,L);
    ksi=zeros(N,1);
    z=zeros(N,1);

    if L==1
        for k=N:M
            x=signal(k-(0:N-1));
            e(k)=d(k)-x'*w;
            ksi=e(k)*x;
            z=x'*ksi*x;
            mi=ksi'*z/norm(z,2)^2;
            w=w+mi*ksi;
            w_mk_enlms(k)=w_mk_enlms(k)+norm((h-w),2).^2;
        end
    elseif L>=1
        for k=N:M
            x=signal(k-(0:N-1));
            X=[x,X(:,1:L-1)];

            j=0;
            for i=k-L+1:k
                e(i)=d(i)-X(:,end-j)'*w;
                j=j+1;
            end

            j=0;
            for i=k-L+1:k
                ksi=ksi+e(i)*X(:,end-j);
                j=j+1;
            end
            ksi=ksi/L;

            j=0;
            for i=k-L+1:k
                z=z+X(:,end-j)'*ksi*X(:,end-j);
                j=j+1;
            end
            z=z/L;

            mi=ksi'*z/norm(z,2)^2;
            w=w+mi*ksi;
            w_mk_enlms(k)=w_mk_enlms(k)+norm((h-w),2).^2;

        end
    end
    w_mk_enlms=w_mk_enlms/length(w_mk_enlms);

    %LMS
    mi = 0.01;
    H = zeros(N,1);
    e1 = zeros(M,1);

    for n = N:M
        X = signal(n:-1:n-(N-1));
        y = H'*X;
        e1(n) = d(n)-y;
        H = H + mi*e1(n)*X;
        w_mk_lms(n)=w_mk_lms(n)+norm((h-H),2).^2;
    end
     w_mk_lms=w_mk_lms/length(w_mk_lms);


    %RLS
    H = zeros(N,1);
    lambda = 0.9984;

    e2 = zeros(M,1);
    C=3.2*eye(N);

    for n = 1:M
        if n>=N
            X = signal(n+(0:-1:-N+1));
        else
            X = [signal(n:-1:1);zeros(N-n,1)];
        end
        y = H'*X;
        e2(n) = d(n)-y;
        mi=X'*C*X;
        g = C*X/(lambda+mi);
        C=1/lambda*C-1/lambda*g*X'*C;
        H = H+g*e2(n);
        w_mk_rls(n)=w_mk_rls(n)+norm((h-H),2).^2;

    end
    w_mk_rls=w_mk_rls/length(w_mk_rls);

    %NLMS
    alfa=1e-5;
    mi_norm=1.45;
    H=zeros(N,1);
    e3=zeros(M,1);

    for n = N:M
        X=signal(n:-1:n-(N-1));
        y=H'*X;
        e3(n)=d(n)-y;
        mi=mi_norm/(alfa+X'*X);
        H=H+mi*e3(n)*X;
        w_mk_nlms(n)=w_mk_nlms(n)+norm((h-H),2).^2;

    end
    w_mk_nlms=w_mk_nlms/length(w_mk_nlms);


end

MSD_enlms = w_mk_enlms/iter;
MSD_lms = w_mk_lms/iter;
MSD_rls = w_mk_rls/iter;
MSD_nlms = w_mk_nlms/iter;


figure(1);
semilogy(MSD_enlms.^2)
hold on
semilogy(MSD_lms.^2)
hold on
semilogy(MSD_rls.^2)
hold on
semilogy(MSD_nlms.^2)
title('Mean-Square Deviation (MSD) in dB');
xlabel('Iteration');
ylabel('MSD (dB)');
legend('ENLMS', 'LMS', 'RLS', 'NLMS')