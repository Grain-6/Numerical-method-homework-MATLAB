clc
clear all
syms x;

f = 2 * cos(2*x) + 4 * x * sin(x) - x^2 - 2; %funtion
p0 = pi/2; 
tol = 1E-5; % accuracy
maxK = 1E3; % maximum iterations allowed

[pn,kn,Yn]=NTM(f,p0,tol,maxK);

fprintf('\n')

f = 2 * cos(2*x) + 4 * x * sin(x) - x^2 - 2; %funtion
p0 = pi/2;  p1 = 1;
tol = 1E-5; % accuracy
maxK = 1E3; % maximum iterations allowed

[ps,ks,Ys]=SEM(f,p0,p1,tol,maxK);


function [p,k,Y]=NTM(f,p0,tol,maxK)
%p0:starting point
%f:equation to solve(function=0)
%maxK:maximum iterations allowed
%tolr:tolerance
%k:times of iteration
%p:result
    syms x;
    P(1)=p0;
    k=2;
    df=diff(f);     %the f'(x) calculated by diff()
   
    P(k)=P(k-1)-subs(f,x,P(k-1))/subs(df,x,P(k-1)); %second iteration
    
    while k<=maxK
        err=abs(P(k)-P(k-1));    %error value
        if(err<tol)
            fprintf("After %d iterations, result reached accurancy of %f by Newton's method.\n",k-1,tol)
            break;
        end
        k=k+1;
        P(k)=P(k-1)-subs(f,x,P(k-1))/subs(df,x,P(k-1));  %iterate by Newton's method 
    end      %k-1 iterations in total
    if(k-1==maxK) 
        disp("Exceeded the maximum number of iterations!");
    end
    p=P(k); 
    k=k-2;
    Y=P;
    fprintf('Starting point = %f \n',p0);
    fprintf('Solution = %f \n',p);
end


function [p,k,Y]=SEM(f,p0,p1,tol,maxK)
%p0,p1:starting point
%f:equation to solve(function=0)
%maxK:maximum iterations allowed
%tolr:tolerance
%k:times of iteration
%p:result

    syms x;
    P(1)=p0;
    P(2)=p1;
    k=3;
    
    df=diff(f);     %利the f'(x) calculated by diff()
    
    P(k)=P(k-1)-subs(f,x,P(k-1))*((P(k-1) - P(k-2))/(subs(df,x,P(k-1)-subs(df,x,P(k-2))))); %second iteration
    
    while k<=maxK
        err=abs(P(k)-P(k-1));    %error value
        if(err<tol)
            fprintf("After %d iterations, result reached accurancy of %f by the secant method.\n",k-2,tol)
            break;
        end
        k=k+1;
        P(k)=P(k-1)-subs(f,x,P(k-1))/subs(df,x,P(k-1));   %iterate by the secant method 
    end      %k-2 iterations in total
    if(k-2==maxK) 
        disp("Exceeded the maximum number of iterations!");
    end
    p=P(k); 
    k=k-2;
    Y=P;
    fprintf('Starting point = %f,%f \n',p0,p1);
    fprintf('Solution = %f \n',p);
end
