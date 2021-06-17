% La regla de Simpson comienza construyendo un polinomio de segundo grado 
% que pasa por los puntos (xi,f(xi)),(xi+h/2,f(xi+h/2)) y (xi+h,f(xy0+h)) 
% = (xy0+1,f(xy0+1)) y usando el área debajo de ese polinomio para estimar 
% el área de la curva. Los cálculos algebraicos involucrados en la 
% construcción del polinomio y el cálculo del área bajo el polinomio son 
% algo complicados, por lo que usaré el paquete sympy para realizar los 
% cálculos de álgebra necesarios.
% La suma de todas estas estimaciones de área sobre todos los subintervalos 
% produce una estimación para el área bajo la curva


f = @(x) (sqrt(1+cos(x.^2)));%funcion a calcular 
a = 1;%intervalo inferior de la integral
b = pi;%intervalo superior de la integral
tol = 0.0001;%error absoluto permitido
S0 = 0;
N = 20;%intervalos
j = 0;
h = (b-a)/2;
S = f(h)*(b-a);
while (abs(S-S0)>tol)
S0 = S;
h = (b-a)/(2*N);
i = 0:N-1;
xi = a+2*i*h;
xi1 = a+2*(i+0.5)*h;
xi2 = a+2*(i+1)*h;
S = (h/3)*sum(f(xi)+4*f(xi1)+f(xi2));%metodo simpson
j = j + 1;
N = 2.^j;
end
S