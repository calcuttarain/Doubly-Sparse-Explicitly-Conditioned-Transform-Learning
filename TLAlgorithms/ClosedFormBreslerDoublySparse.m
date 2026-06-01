function [B,XB,error]= ClosedFormBreslerDoublySparse(B, Y, Y_test, numiter, l2, l3, T1, STY, STY_te);

cbb = 0;

%This is an implementation of the transform learning algorithm with closed-form solutions for the sparse coding and transform update steps that was presented in the following papers:
%1) S. Ravishankar and Y. Bresler, "Closed-form solutions within sparsifying transform learning," in Proc. IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 2013, pp. 5378-5382.
%2) S. Ravishankar and Y. Bresler, Closed-Form Optimal Updates In Transform Learning, in Signal Processing with Adaptive Sparse Structured Representations (SPARS) workshop, July 2013.
%3) S. Ravishankar and Y. Bresler, \ell_0 Sparsifying transform learning with efficient optimal updates and convergence guarantees, IEEE Transactions on Signal Processing, vol. 63, no. 9, pp. 2389-2404, May 2015.

%We employ alternating minimization here to solve a transform learning problem that involves a constraint on the adaptive transform domain sparsity of each training signal, 
%and a transform learning regularizer that involves a log-determinant and a Frobenius norm penalty.
%The algorithm iterates over a sparse coding step and a transform update step, both of which involve efficient update procedures.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inputs: 1) W : Initial Transform
%        2) Y : Training Matrix with signals as columns
%        3) numiter:  Number of iterations of alternating minimization
%        4) l2 : Weight for log-determinant penalty
%        5) l3 : Weight for Frobenius norm penalty
%        6) STY: Vector containing maximum allowed sparsity levels for each training signal.

%Outputs:  1) W: Learnt Transform
%          2) XB: Learnt Sparse Code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial steps
[K,n]=size(B);XB=zeros(K,size(Y,2));
[U,S,V]=svd((Y*Y') + (l3*eye(n)));
LL=U*(S^(1/2))*V';
LL2=(inv(LL));


ix=find(STY>0); q=Y(:,ix); STY=STY(:,ix); N=size(q,2);
ez=K*(0:(N-1));STY=STY + ez;
Y = Y(:, ix);

ix_te=find(STY_te>0); q_te=Y_test(:,ix_te); STY_te=STY_te(:,ix_te); N_te=size(q_te,2);
ez_te=K*(0:(N_te-1));STY_te=STY_te + ez_te;
Y_test = Y_test(:, ix_te);

error.m1.tr = []; error.m1.te = [];
error.m2.tr = []; error.m2.te = [];

%Algorithm iterations in a FOR Loop
for i=1:numiter
    
    %Sparse Coding Step
    X1=B*q;
    [s]=sort(abs(X1),'descend');
    X = X1.*(bsxfun(@ge,abs(X1),s(STY)));

    X1_te=B*q_te;
    [s_te]=sort(abs(X1_te),'descend');
    X_test = X1_te.*(bsxfun(@ge,abs(X1_te),s_te(STY_te)));
%     X = X/norm(X, 'fro');
    
    %Transform Update Step
    [Q1,Si,R]=svd(LL2*q*X');
    sig=diag(Si);
    gamm=(1/2)*(sig + (sqrt((sig.^2) + 2*l2)));
    D=R*(diag(gamm))*Q1';
    B=D*(LL2);
%     W = normc(W);
%     W = sqrt(n)*W/norm(W, 'fro');

    %Post-thresholding
    if(i>cbb) %done only after initial `cbb' number of iterations
        [~,vu]=sort(abs(B(:)),'descend');
        Bg=zeros(size(B));
        Bg(vu(1:T1))=B(vu(1:T1));
        B=Bg;
    end

    %check for zero determinant
    % cll=1e-1;
    % while(abs(det(B)) <= 10^(-250))
    %     B = B + ((rand(K,n) - 0.5)*cll);
    % end

    error.m1.tr = [error.m1.tr norm(X - B*Y, 'fro')];
    error.m2.tr = [error.m2.tr norm(X - B*Y, 'fro')/norm(B*Y, 'fro')];

    error.m1.te = [error.m1.te norm(X_test - B*Y_test, 'fro')];
    error.m2.te = [error.m2.te norm(X_test - B*Y_test, 'fro')/norm(B*Y_test, 'fro')];
end
XB(:,ix)=X;
