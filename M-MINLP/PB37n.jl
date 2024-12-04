# Import necessary libraries
import DataFrames
import LinearAlgebra
# Create system impedance matrices
function Matriz36Z(C) #Ω/km
    if C == 1
        Zm = [0.2926 + im*0.1973   0.0673 - im*0.0368   0.0337 - im*0.0417;
              0.0673 - im*0.0368   0.2646 + im*0.1900   0.0673 - im*0.0368;
              0.0337 - im*0.0417   0.0673 - im*0.0368   0.2926 + im*0.1973];
    elseif C == 2
        Zm = [0.4751 + im*0.2973   0.1629 - im*0.0326   0.1234 - im*0.0607;
              0.1629 - im*0.0326   0.4488 + im*0.2678   0.1629 - im*0.0326;
              0.1234 - im*0.0607   0.1629 - im*0.0326   0.4751 + im*0.2973];
    elseif C == 3
        Zm = [1.2936 + im*0.6713   0.4871 + im*0.2111   0.4585 + im*0.1521;
              0.4871 + im*0.2111   1.3022 + im*0.6326   0.4871 + im*0.2111
              0.4585 + im*0.1521   0.4871 + im*0.2111   1.2936 + im*0.6713];
    elseif C == 4
        Zm =  [2.0952 + im*0.7758   0.5204 + im*0.2738   0.4926 + im*0.2123;
               0.5204 + im*0.2738   2.1068 + im*0.7398   0.5204 + im*0.2738;
               0.4926 + im*0.2123   0.5204 + im*0.2738   2.0952 + im*0.7758];
    end
end
# Define base values 
Vb = 4.8/sqrt(3);  # kV
Sb = 1000; # kVA
Ib = Sb/Vb; # A
Zb = ((1000*Vb)^2)/(Sb*1000); # Ω
# Branch information
lineas = DataFrames.DataFrame([
    (1,2,1,1850),
    (2,3,2,960),
    (3,24,4,400),
    (3,27,3,360),
    (3,4,2,1320),
    (4,5,4,240),
    (4,9,3,600),
    (5,6,3,280),
    (6,7,4,200),
    (6,8,4,280),
    (9,10,3,200),
    (10,23,3,600),
    (10,11,3,320),
    (11,13,3,320),
    (11,12,4,320),
    (13,14,3,560),
    (14,18,3,640),
    (14,15,4,520),
    (15,16,4,200),
    (15,17,4,1280),
    (18,19,3,400),
    (19,20,3,400),
    (20,22,3,400),
    (20,21,4,200),
    (24,26,4,320),
    (24,25,4,240),
    (27,28,3,520),
    (28,29,4,80),
    (28,31,3,800),
    (29,30,4,520),
    (31,34,4,920),
    (31,32,3,600),
    (32,33,4,280),
    (34,36,4,760),
    (34,35,4,120)
]); 
DataFrames.rename!(lineas, [:i, :j, :Zij, :Lij]);
# Node information
nodos = DataFrames.DataFrame([
    (1,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (2,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),140,70,140,70,350,175,0),
    (3,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (4,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (5,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,42,21,0),
    (6,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),42,21,0,0,0,0,0),
    (7,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),42,21,42,21,42,21,0),
    (8,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),42,21,0,0,0,0,0),
    (9,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,85,40,0),
    (10,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (11,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (12,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,42,21,0),
    (13,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),85,40,0,0,0,0,0),
    (14,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,42,21,0),
    (15,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (16,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,85,40,0),
    (17,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,42,21,0,0,0),
    (18,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),140,70,0,0,0,0,0),
    (19,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),126,62,0,0,0,0,0),
    (20,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (21,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,85,40,0),
    (22,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,42,21,0),
    (23,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,85,40,0,0,0),
    (24,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (25,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,85,40,0),
    (26,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),8,4,85,40,0,0,0),
    (27,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,85,40,0),
    (28,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (29,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),17,8,21,10,0,0,0),
    (30,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),85,40,0,0,0,0,0),
    (31,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,85,40,0),
    (32,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (33,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,42,21,0,0,0),
    (34,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,0,0,0,0,0),
    (35,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,140,70,21,10,0),
    (36,1,exp(im*deg2rad(-120)),exp(im*deg2rad(120)),0,0,42,21,0,0,0),
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
    local Yp3[3*k-2:3*k,3*k-2:3*k] = inv((Matriz36Z(lineas.Zij[k])*lineas.Lij[k]*0.000189394)/Zb);
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
#3. Define restrictions  
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
@objective(OPF,Min,real(transpose(V)*conj(Ybus3*V)));
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



