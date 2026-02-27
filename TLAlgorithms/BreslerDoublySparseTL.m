function [B,XB,i,error,error2]= BreslerDoublySparseTL(B,Y,numiter,l2,l3,T1,mu,numgrad,STY,cbb,stopcn,stopth)

%This is an implementation of a transform learning algorithm that was presented in the following papers:
%1) S. Ravishankar and Y. Bresler, Learning doubly sparse transforms for images, IEEE Transactions on Image Processing, vol. 22, no. 12, pp. 4598-4612, Dec. 2013.
%2) S. Ravishankar and Y. Bresler, Learning doubly sparse transforms for image representation, in IEEE International Conference on Image Processing (ICIP), 2012, pp. 685-688.

%We employ alternating minimization here to solve a transform learning problem that involves a constraint on the adaptive transform domain sparsity of each training signal, 
%a constraint on the sparsity of the matrix B in  the decomposition W = B*\Phi, and a transform learning regularizer that involves a log-determinant penalty and a Frobenius norm penalty.
%The algorithm iterates over a sparse coding step and a transform update step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inputs: 1) B : Initialization with T1 non-zeros. Note that B is the sparse Matrix in the decomposition W=B*\Phi
%        2) Y : Training Matrix (with signals as columns) that has been pre-processed by applying the fixed transform \Phi to it.
%        3) numiter:  Number of iterations of alternating minimization
%        4) l2 : Weight for the log-determinant penalty in the problem formulation
%        5) l3 : Weight for the Frobenius norm penalty in the problem formulation
%        6) T1  : maximum number of non-zeros allowed in the sparse matrix B
%        7) mu: Step size within CG in the transform update step
%        8) numgrad: Number of iterations of CG within the transform update step
%        9) STY: Vector containing maximum allowed sparsity levels for each training signal
%       10) cbb: Number of initial algorithm iterations with no post-thresholding of the transform in the transform update step
%       11) stopcn: If set to 1, the relative iterate change is computed and monitored over the algorithm iterations
%       12) stopth: Threshold for the relative iterate change. Algorithm iterations terminate when the relative iterate change drops below `stopth'.
%                   This parameter is only used by the algorithm when `stopcn' is set to 1.

%Outputs:  1) B: Learnt Sparse matrix B in the decomposition W=B*\Phi
%          2) XB: Matrix with learnt sparse representations of training data as columns
%          3) i: Iteration number for the last executed algorithm iteration.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial steps
[K,n]=size(B);XB=zeros(K,size(Y,2)); Bgh=B;

ix=find(STY>0); q=Y(:,ix); STY=STY(:,ix); N=size(q,2); 
ez=K*(0:(N-1));STY=STY + ez;

ZX=Y*Y'; %pre-computed for transform update step

error = [];
error2 = [];

%Algorithm iterations in a FOR Loop
for i=1:numiter
 
    %%%%%%Sparse Coding Step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X1=(sparse(B))*q;
    [s]=sort(abs(X1),'descend');
    X = X1.*(bsxfun(@ge,abs(X1),s(STY)));
           
    %%%%%%Transform Update Step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           RSu=(q*X')';
           
           %CG iterations
           for j=1:numgrad

               cc = (-2*RSu) + (2*B*ZX');
               cc2= (-1)*((inv(B))');
               cc3=2*B;              
               deD= l2*cc2 + cc + l3*cc3;  %gradient of objective
                                            
               if(j==1)
                   g=deD;
                   d=-g;
               else
                   be = ((norm(deD,'fro'))^2)/((norm(g,'fro'))^2);
                   g=deD;
                   d= - g + be*d;
               end
    
               B= B + (mu)*d;   %CG with constant step size

           end
           
           %Post-thresholding of transform
           if(i>cbb) %done only after initial `cbb' number of iterations
           [~,vu]=sort(abs(B(:)),'descend');
           Bg=zeros(size(B));
           Bg(vu(1:T1))=B(vu(1:T1));
           B=Bg;
           end

           %check for zero determinant
           cll=1e-1;  
           while(abs(det(B)) <= 10^(-250))
             B = B + ((rand(K,n) - 0.5)*cll);               
           end

           error = [error norm(X - B*q, 'fro')];
           error2 = [error2 norm(X - B*q, 'fro')/norm(B*q, 'fro')];
          
           
           
%Monitor the relative iterate change when the stopping criterion parameter is active.
if(stopcn==1)
   sh=(norm(B-Bgh,'fro'))/(norm(Bgh,'fro'));
   Bgh=B;
end

%Check the stopping condition when the stopping criterion parameter is active.
if(stopcn==1)
if(i>cbb)
if(sh<stopth)
    break;  %break from loop if stopping threshold satisfied.
end
end
end
    
end
XB(:,ix)=X;

