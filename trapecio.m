%La integral definida de f(x) es igual al area neta debajo de la curva y
%=f(x) durante el intervalo [a,b]. Riemman suma integrales definidas aproximadas usando sumas
% de rectangulos para aproximar el area La regla del trapezoide da una mejor aproximación 
% de una integral definida sumando las áreas de los trapezoides que conectan los puntos.
% (xi-1,0),(xi,0)(xi-1,f(xi-1)),(xi,f(xi))
% para cada subintervalo [xi-1,xi] de una particion, tenga en cuenta que el
% area de cada trapezoide es la suma de un rectangulo y un triangulo 
% (xi-xi-1)f(xi-1)+(1/2)(xi-xi-1)(f(xi)-f(xi-1)),(xi,f(xi))
function R = trapecio(f,a,b,depth,tol)
%f es la funcion objetivo
%a & b son los limites de integracion
%depth es el numero de intervalos
%tol es el error minimmo permitido en este caso 0.001 que es menor a 10^-2
M=1;
h=b-a;
err=1;
J=0;
R=zeros(depth,depth);
%Inicio de la regla trapezoidal compuesta
% y evaluación de los valores de la primera columna
R(1,1)=h*(feval(f,a)+feval(f,b))/2;
while((err>tol)&&(J<depth))||(J<depth)
    J=J+1;
    h=h/2;
    sum=0;
    for p=1:M
        x=a+h*(2*p-1);
        sum=sum+feval(f,x);
    end
    R(J+1,1)=R(J,1)/2+h*sum;
    M=2*M;
%fin de la regla trapezoidal
% y evaluacion de las otras columnas
    for K=1:J
        R(J+1,K+1) = R(J+1,K)+(R(J+1,K)-R(J,K))/(4^K-1);
    end
    err = abs(R(J,J)-R(J+1,K+1));
end
%