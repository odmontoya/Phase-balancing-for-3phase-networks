# Import necessary libraries
import DataFrames
import LinearAlgebra
# Create system impedance matrices
function Matriz25Z(C) #Ω/km
    if C == 1
        Zm = [0.3686 + im*0.6852   0.0169 + im*0.1515   0.0155 + im*0.1098;
              0.0169 + im*0.1515   0.3757 + im*0.6715   0.0188 + im*0.2072;
              0.0155 + im*0.1098   0.0188 + im*0.2072   0.3723 + im*0.6782];
    elseif C == 2
        Zm = [0.9775 + im*0.8717   0.0167 + im*0.1697   0.0152 + im*0.1264;
              0.0167 + im*0.1697   0.9844 + im*0.8654   0.0186 + im*0.2275;
              0.0152 + im*0.1264   0.0186 + im*0.2275   0.981 + im*0.8648];
    elseif C == 3
        Zm = [1.928 + im*1.4194    0.0161 + im*0.1183   0.0161 + im*0.1183;
              0.0161 + im*0.1183   1.9308 + im*1.4215   0.0161 + im*0.1183;
              0.0161 + im*0.1183   0.0161 + im*0.1183   1.9337 + im*1.4236];
    end
end
# Define base values 
Vb = 4.16/sqrt(3); # kV
Sb = 1000; # kVA
Ib = Sb/Vb; # A
Zb = ((1000*Vb)^2)/(Sb*1000); # Ω
# Branch information
lineas = DataFrames.DataFrame([
    (1,2,1,0.189394),
    (2,3,1,0.094697),
    (2,6,2,0.094697),
    (3,4,1,0.094697),
    (3,18,2,0.094697),
    (4,5,2,0.094697),
    (4,23,2,0.0757576),
    (6,7,2,0.094697),
    (6,8,2,0.189394),
    (7,9,2,0.094697),
    (7,14,2,0.094697),
    (7,16,2,0.094697),
    (9,10,2,0.094697),
    (10,11,2,0.0568182),
    (11,12,3,0.0378788),
    (11,13,3,0.0378788),
    (14,15,2,0.0568182),
    (14,17,3,0.0568182),
    (18,20,2,0.094697),
    (18,21,3,0.0757576),
    (20,19,3,0.0757576),
    (21,22,3,0.0757576),
    (23,24,2,0.0757576),
    (24,25,3,0.0757576),
]); 
DataFrames.rename!(lineas, [:i, :j, :Zij, :Lij]);
# Node information
nodos = DataFrames.DataFrame([
    (1,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (2,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (3,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),36,21.6,28.8,19.2,42,26.4,0),
    (4,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),57.6,43.2,4.8,3.4,48,30,0),
    (5,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),43.2,28.8,28.8,19.2,36,24,0),
    (6,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),43.2,28.8,33.6,24,30,30,0),
    (7,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (8,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),43.2,28.8,28.8,19.2,3.6,2.4,0),
    (9,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),72,50.4,38.4,28.8,48,30,0),
    (10,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),36,21.6,28.8,19.2,32,26.4,0),
    (11,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),50.4,31.7,24,14.4,36,24,0),
    (12,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),57.6,36,48,33.6,48,36,0),
    (13,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),64.8,21.6,33.6,21.1,36,24,0),
    (14,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),57.6,36,38.4,28.8,60,42,0),
    (15,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),7.2,4.3,4.8,2.9,6,3.6,0),
    (16,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),57.6,4.3,3.8,28.8,48,36,0),
    (17,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),57.6,43.2,33.6,24,54,38.4,0),
    (18,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),57.6,43.2,38.4,28.8,48,36,0),
    (19,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),8.6,6.5,4.8,3.4,6,4.8,0),
    (20,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),50.4,36,38.4,28.8,54,38.4,0),
    (21,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),5.8,4.3,3.4,2.4,5.4,3.8,0),
    (22,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),72,50.4,57.6,43.2,60,48,0),
    (23,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),8.6,64.8,4.8,3.8,60,42,0),
    (24,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),50.4,36,43.2,30.7,4.8,3.6,0),
    (25,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),8.6,6.5,4.8,2.9,6,4.2,0),
    
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
    local Yp3[3*k-2:3*k,3*k-2:3*k] = inv((Matriz25Z(lineas.Zij[k])*lineas.Lij[k])/Zb);
end
Ybus3 = A3*Yp3*transpose(A3);
# Model information
M = [1 -1 0; 0 1 -1; -1 0 1];
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
set_attribute(OPF, "honor_original_bounds", "yes")
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
@variable(OPF, X[k in 1:3*NN, j in 1:3], Bin);
#3. Define constraints  
@constraint(OPF, V[3*slack-2] == 1.0 + im*0.0);
@constraint(OPF, V[3*slack-1] == -0.5 - im*0.866025403784439);
@constraint(OPF, V[3*slack] == -0.5 + im*0.866025403784439);
@constraint(OPF, X[3*slack-2:3*slack,:] .== [1 0 0; 0 1 0; 0 0 1]);
for k = 1:NN
    @constraint(OPF, Ig[3*k-2:3*k] - Idy[3*k-2:3*k] - transpose(M)*Idd[3*k-2:3*k] == (Ybus3[3*k-2:3*k,:]*V));
    @constraint(OPF, conj(Sg[3*k-2:3*k]) == LinearAlgebra.diagm(conj(V[3*k-2:3*k]))*Ig[3*k-2:3*k]);
    @constraint(OPF, X[3*k-2:3*k,:]*conj(Sdy[3*k-2:3*k]) == LinearAlgebra.diagm(conj(V[3*k-2:3*k]))*Idy[3*k-2:3*k]);
    @constraint(OPF, X[3*k-2:3*k,:]*conj(Sdd[3*k-2:3*k]) == LinearAlgebra.diagm(conj(M*V[3*k-2:3*k]))*Idd[3*k-2:3*k]);
    @constraint(OPF, sum(X[3*k-2,:]) == 1);
    @constraint(OPF, sum(X[3*k-1,:]) == 1);
    @constraint(OPF, sum(X[3*k,:]) == 1);
    @constraint(OPF, sum(X[3*k-2:3*k,1]) == 1);
    @constraint(OPF, sum(X[3*k-2:3*k,2]) == 1);
    @constraint(OPF, sum(X[3*k-2:3*k,3]) == 1);
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
#ABC
M1 = [1 0 0;  0 1 0; 0 0 1];
#BCA
M2 = [0 1 0;  0 0 1; 1 0 0];
#CAB
M3 = [0 0 1;  1 0 0; 0 1 0];
#ACB
M4 = [1 0 0;  0 0 1; 0 1 0];
#CBA
M5 = [0 0 1;  0 1 0; 1 0 0];
#BAC
M6 = [0 1 0;  1 0 0; 0 0 1];
sol = zeros(1,NN)
for k = 1:NN
    if round.(value.(X[3*k-2:3*k,:])) == M1
        sol[1,k] = 1;
    elseif round.(value.(X[3*k-2:3*k,:])) == M2
        sol[1,k] = 2;
    elseif round.(value.(X[3*k-2:3*k,:])) == M3
        sol[1,k] = 3;  
    elseif round.(value.(X[3*k-2:3*k,:])) == M4
        sol[1,k] = 4;
    elseif round.(value.(X[3*k-2:3*k,:])) == M5
        sol[1,k] = 5;
    elseif round.(value.(X[3*k-2:3*k,:])) == M6
        sol[1,k] = 6;
    end  
end


