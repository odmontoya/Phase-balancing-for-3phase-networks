# Import necessary libraries
import DataFrames
import LinearAlgebra
# Create system impedance matrices
function Matriz8Z(C) #Ω/km
    if C == 1
        Zm = [0.093654 + im*0.040293    0.031218 + im*0.013431	0.031218 + im*0.013431;
              0.031218 + im*0.013431	0.093654 + im*0.040293	0.031218 + im*0.013431;
              0.031218 + im*0.013431	0.031218 + im*0.013431	0.093654 + im*0.040293];
    elseif C == 2
        Zm = [0.15609 + im*0.067155	0.05203 + im*0.022385	0.05203 + im*0.022385;
              0.05203 + im*0.022385	0.15609 + im*0.067155	0.05203 + im*0.022385;
              0.05203 + im*0.022385	0.05203 + im*0.022385	0.15609 + im*0.067155];
    elseif C == 3
        Zm = [0.046827 + im*0.0201465	0.015609 + im*0.0067155	 0.015609 + im*0.0067155;
              0.015609 + im*0.0067155	0.046827 + im*0.0201465	 0.015609 + im*0.0067155;
              0.015609 + im*0.0067155	0.015609 + im*0.0067155	 0.046827 + im*0.0201465];
    elseif C == 4
        Zm = [0.031218 + im*0.013431	0.010406 + im*0.004477	0.010406 + im*0.004477;
              0.010406 + im*0.004477	0.031218 + im*0.013431	0.010406 + im*0.004477;
              0.010406 + im*0.004477	0.010406 + im*0.004477	0.031218 + im*0.013431];
    elseif C == 5
        Zm = [0.062436 + im*0.026862	0.020812 + im*0.008954	0.020812 + im*0.008954;
              0.020812 + im*0.008954	0.062436 + im*0.026862	0.020812 + im*0.008954;
              0.020812 + im*0.008954	0.020812 + im*0.008954	0.062436 + im*0.026862];
    elseif C == 6
        Zm = [0.078045 + im*0.0335775	0.026015 + im*0.0111925	 0.026015 + im*0.0111925;
              0.026015 + im*0.0111925	0.078045 + im*0.0335775	 0.026015 + im*0.0111925;
              0.026015 + im*0.0111925	0.026015 + im*0.0111925	 0.078045 + im*0.0335775];
    end
end
# Define base values
Vb = 11/sqrt(3); # kV
Sb = 1000; # kVA
Ib = Sb/Vb; # A
Zb = ((1000*Vb)^2)/(Sb*1000); # Ω
# Branch information
lineas = DataFrames.DataFrame([
    (1,	2,	1,	1),
    (2,	3,	2,	1),
    (2,	5,	3,	1),
    (2,	7,	3,	1),
    (3,	4,	4,	1),
    (3,	8,	5,	1),
    (5,	6,	6,	1),
]); 
DataFrames.rename!(lineas, [:i, :j, :Zij, :Lij]);
# Node information
nodos = DataFrames.DataFrame([
    (1,   1,   exp(im*deg2rad(-120)), exp(im*deg2rad(120)), 0,    0,    0,     0,     0,     0,     0),
    (2,   1,   exp(im*deg2rad(-120)), exp(im*deg2rad(120)), 519,  250,	259,   126,	  515,	 250,   0),
    (3,   1,   exp(im*deg2rad(-120)), exp(im*deg2rad(120)), 0,	  0,	259,   126,	  486,	 235,   0),
    (4,   1,   exp(im*deg2rad(-120)), exp(im*deg2rad(120)), 0,	  0,	0,	   0, 	  324,	 157,   0),
    (5,   1,   exp(im*deg2rad(-120)), exp(im*deg2rad(120)), 0,	  0,	0,	   0,	  226,	 109,   0),
    (6,   1,   exp(im*deg2rad(-120)), exp(im*deg2rad(120)), 0,	  0,	0,	   0,	  145,	 70,    0),
    (7,   1,   exp(im*deg2rad(-120)), exp(im*deg2rad(120)), 486,  235,	0,	   0,	  0,	 0,     0),
    (8,   1,   exp(im*deg2rad(-120)), exp(im*deg2rad(120)), 0,	  0,    267,   129,	  0,	 0,     0),
]);
DataFrames.rename!(nodos, [:i, :Vai0, :Vbi0, :Vci0, :Pai, :Qai, :Pbi, :Qbi, :Pci, :Qci, :Type])
#Type = 0 ---> Y
#Type = 1 ---> Δ
# Convert to per unit
nodos.Pai = nodos.Pai/Sb;
nodos.Qai = nodos.Qai/Sb;
nodos.Pbi = nodos.Pbi/Sb;
nodos.Qbi = nodos.Qbi/Sb;
nodos.Pci = nodos.Pci/Sb;
nodos.Qci = nodos.Qci/Sb;
# Ybus matrix formation
NN = size(nodos,1);
NL = size(lineas,1);
A3 = zeros(3*NN,3*NL);
Yp3 = complex(zeros(3*NL,3*NL));
for k = 1:NL
    Ni = lineas.i[k];
    Nj = lineas.j[k];
    A3[3*Ni-2:3*Ni,3*k-2:3*k] = [1 0 0; 0 1 0; 0 0 1];
    A3[3*Nj-2:3*Nj,3*k-2:3*k] = [-1 0 0; 0 -1 0; 0 0 -1];
    local Yp3[3*k-2:3*k,3*k-2:3*k] = inv((Matriz8Z(lineas.Zij[k])*lineas.Lij[k])/Zb);
end
Ybus3 = A3*Yp3*transpose(A3);
# Model information
# Matrix M
M = [1 -1 0; 0 1 -1; -1 0 1];
#Conection matrices
#ABC
M1 = [1 + 0im 0 0;  0 1 0; 0 0 1];
#BCA
M2 = [0 + 0im 1 0;  0 0 1; 1 0 0];
#CAB
M3 = [0 + 0im 0 1;  1 0 0; 0 1 0];
#ACB
M4 = [1 + 0im 0 0;  0 0 1; 0 1 0];
#CBA
M5 = [0 + 0im 0 1;  0 1 0; 1 0 0];
#BAC
M6 = [0 + 0im 1 0;  1 0 0; 0 0 1];
H = [M1;M2;M3;M4;M5;M6];
slack = 1;
Vmin = 0.90;
Vmax = 1.10;
# Wye loads
Sdy = complex(zeros(3*(NN),1));
# Triangle loads
Sdd = complex(zeros(3*(NN),1));
for k = 1:NN
    if nodos.Type[k] == 0
        Sdy[3*k-2:3*k,1] = [nodos.Pai[k] + im*nodos.Qai[k];
                            nodos.Pbi[k] + im*nodos.Qbi[k];
                            nodos.Pci[k] + im*nodos.Qci[k]]; 
    else    
        Sdd[3*k-2:3*k,1] = [nodos.Pai[k] + im*nodos.Qai[k];
                            nodos.Pbi[k] + im*nodos.Qbi[k];
                            nodos.Pci[k] + im*nodos.Qci[k]]; 
    end
end
# Optimization model
# 1. Create an optimization model called OPF using the Bonmin solver
using JuMP
import AmplNLWriter
import Bonmin_jll
OPF = Model(() -> AmplNLWriter.Optimizer(Bonmin_jll.amplexe))
set_attribute(OPF, "bonmin.nlp_log_level", 0)
# 2. Define optimization variables 
@variable(OPF,V[k in 1:3*NN] in ComplexPlane());
for k = 1:NN 
    set_start_value(real(V[3*k-2]),1.0);
    set_start_value(real(V[3*k-1]),-0.5);
    set_start_value(real(V[3*k]),-0.5);
    set_start_value(imag(V[3*k-2]),0.0);
    set_start_value(imag(V[3*k-1]),-0.866025403784439);
    set_start_value(imag(V[3*k]),0.866025403784439);
end
@variable(OPF, Sg[k in 1:3*NN] in ComplexPlane());
@variable(OPF, Ig[k in 1:3*NN] in ComplexPlane());
@variable(OPF, Idy[k in 1:3*NN] in ComplexPlane());
@variable(OPF, Idd[k in 1:3*NN] in ComplexPlane());
@variable(OPF, X[k in 1:NN, j in 1:6], Bin);
#3. Define constraints
@constraint(OPF, V[3*slack-2] == 1.0 + im*0.0);
@constraint(OPF, V[3*slack-1] == -0.5 - im*0.866025403784439);
@constraint(OPF, V[3*slack] == -0.5 + im*0.866025403784439);
for k = 1:NN
    @constraint(OPF, Ig[3*k-2:3*k] - Idy[3*k-2:3*k] - transpose(M)*Idd[3*k-2:3*k] == (Ybus3[3*k-2:3*k,:]*V));
    @constraint(OPF, conj(Sg[3*k-2:3*k]) == LinearAlgebra.diagm(conj(V[3*k-2:3*k]))*Ig[3*k-2:3*k]);
    @constraint(OPF, sum(X[k,j]*H[3*j-2:3*j,:] for j in 1:6)*conj(Sdy[3*k-2:3*k])  == LinearAlgebra.diagm(conj(V[3*k-2:3*k]))*Idy[3*k-2:3*k]);
    @constraint(OPF, sum(X[k,j]*H[3*j-2:3*j,:] for j in 1:6)*conj(Sdd[3*k-2:3*k]) == LinearAlgebra.diagm(conj(M*V[3*k-2:3*k]))*Idd[3*k-2:3*k]);
    @constraint(OPF, sum(X[k,j] for j in 1:6) == 1);
    for j = 0:2
        @constraint(OPF, abs2(Vmin) <= abs2(V[3*k-j]) <= abs2(Vmax));
    end
    if k != slack 
        @constraint(OPF, Sg[3*k-2:3*k] == 0);
    end

end
# 4. Define the objective function
@objective(OPF,Min,Sb*real(transpose(V)*conj(Ybus3*V)));
# 5. Solve the model
JuMP.optimize!(OPF) 
# 6. Show variables of interest
@show objective_value(OPF);
############### Solution Vector ###################
respuesta = round.(value.(X));
sol = zeros(1,NN);
for k = 1:NN
    a = findall( x -> x == 1, respuesta[k,:]); 
    sol[1,k] = a[1];
end

