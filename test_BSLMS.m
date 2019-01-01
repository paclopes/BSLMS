% test_BSLMS
% 
% Attachment to the Paper
% P. A. C. Lopes, "A Bayesian Step Least Mean Squares Algorithm for
% Gaussian Signals," IEEE Signal Processing Letters, 2019

% simulation parameters
NAvr=50;        % number of Monte Carlo runs
N=16;           % filter size
R=4000;         % simulation time
qv=0.01;        % measurement noise
qu=1;           % reference signal power
sw=1;           % optimal filter coefficient standard deviation
h=zeros(2,1);   % the reference spectral shaping filter
h(1)=1; h(2)=1/2;

% algorithm parameters
alpha = 1/3;            % 1/(eigenvalues spread) should be always < 1
max_trace_A = 1e11;     % maximum value of the trace of A
k = 2;                  % E[e^4]/E[e^2]^2 - 1

% filters
wo=(1.1*sinc((0:N-1)-10.7)-0.7*sinc((0:N-1)-13.2)+0.5*sinc((0:N-1)-5.3))';
wo=sw*wo/sqrt(sum(wo.^2));

% statistics
e_log=zeros(R,NAvr);
d_log=zeros(R,NAvr);
A_log=zeros(R,2,2,NAvr);
b_log=zeros(R,2,NAvr);
qo_log=zeros(R,2,NAvr);
w_msq_log=zeros(R,NAvr);
mu_log=zeros(R,NAvr);

for i=1:NAvr
    rng(217928+i); % seed
    % reference:
    u=sqrt(qu)*conv(randn(1,length(h)+R-1),h)/sqrt(Rm(1,1));u=u(length(h):end);  
    v=sqrt(qv)*randn(R,1);  % measurement noise
    
    % initialization
    w=zeros(N,1); w(1)=1;
    uv=zeros(N,1);
    qo=ones(2,1); % [qw, qv]'
    A=zeros(2);
    b=zeros(2,1);
    
    % simulation
    for n=1:R
        uv=[u(n) uv(1:N-1)']';
        d = uv'*wo + v(n);

        e = d - uv'*w;
        if n>N
            [w, qo, A, b, mu] = BSLMS(uv, d, w, qo, A, b, N, alpha, max_trace_A, k);

            %%%%%% statistics
            e_log(n,i) = e;
            d_log(n,i)= d;
            A_log(n,:,:,i)=A;
            b_log(n,:,i)=b;
            qo_log(n,:,i)=qo;
            w_msq_log(n,i)=mean((w-wo).^2);
            mu_log(n,i)=mu;
        end;
        
    end
end

figure(1);
plot(10*log10(mean(w_msq_log,2)));
xlabel('time');
ylabel('Mean Coefficient Square Error (dBs)');
grid on;

figure(2);
n_t_avr = 10;
desired_av = mean(abs(d_log).^2,2);
desired_av = mean(reshape(desired_av,n_t_avr,[]));
time=n_t_avr/2:n_t_avr:R;
plot(time, 10*log10(desired_av));
hold on;
error_av = mean(abs(e_log).^2,2);
error_av = mean(reshape(error_av,n_t_avr,[]));
plot(time, 10*log10(error_av));
hold off;
legend('desired', 'error');
xlabel('time');
ylabel('Error (dBs)');
grid on;

figure(3);
semilogy(mean(qo_log,3));
hold on;
semilogy(mean(w_msq_log,2));
hold off;
set(gca,'Ylim',[0.000001,1]);
legend('qw','qv', 'true qw');
ylabel('Estimated values');
xlabel('Time');
grid on;

figure(4); 
semilogy(reshape(qo_log(:,2,:),R,[]));
ylabel('Estimates of q_v');
xlabel('Time');
grid on;

figure(5);
semilogy(mean(mu_log,2));
ylabel('\mu');
xlabel('Time');
grid on;

