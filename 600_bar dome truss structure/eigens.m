function [L,X]=eigens(K,M,dof_constrnt)
% PURPOSE
%  Solve the generalized eigenvalue problem
%  [K-LM]X = 0, considering boundary conditions.
%
% INPUT:
%    K : global stiffness matrix, dim(K)= ele_nod x ele_nod
%    M : global mass matrix, dim(M)= ele_nod x ele_nod
%    b : boundary condition matrix
%        dim(dof_active)= dof_active x 1
% OUTPUT:
%    L : eigenvalues stored in a vector with length (ele_nod - dof_active) 
%    X : eigenvectors dim(X)= ele_nod x num_dof, num_dof : number of dof's

  [ele_nod,ele_nod]=size(K);
  dof=(1:ele_nod)';
%
  if nargin==3
    dof_constrnt=dof_constrnt(:);
    dof(dof_constrnt)=[]; 
    if nargout==2
      [X1,D]=eig(K(dof,dof),M(dof,dof));
      [num_dof,num_dof]=size(X1);
      for jj=1:num_dof
        mnorm=sqrt(X1(:,jj)'*M(dof,dof)*X1(:,jj));
        X1(:,jj)=X1(:,jj)/mnorm;
      end
      d=diag(D);
      [L,ii]=sort(d);
      X2=X1(:,ii);
      X=zeros(ele_nod,num_dof);
      X(dof,:)=X2;
    else
      d=eig(K(dof,dof),M(dof,dof));
      L=sort(d);
    end
  else
    if nargout==2
      [X1,D]=eig(K,M);
      for jj=1:ele_nod
        mnorm=sqrt(X1(:,jj)'*M*X1(:,jj));
        X1(:,jj)=X1(:,jj)/mnorm;
      end
      d=diag(D);
      [L,ii]=sort(d);
      X=X1(:,ii);
    else
      d=eig(K,M);
      L=sort(d);
    end
  end

