%%% OwnCode 2
%%% Problem :
%%% min z = 35x1 + 30x2 + 25x3
%%% s.t 
%%%     90x1 + 20x2 + 40x3 >= 200
%%%     30x1 + 80x2 + 60x3 >= 180
%%%     10x1 + 20x2 + 50x3 >= 150
%%%           x1,x2,x3 >= 0
%%% After adding slack variables and artificial variables
%%%     90x1 + 20x2 + 40x3 - s1           + A1           = 200
%%%     30x1 + 80x2 + 60x3      - s2           + A2      = 180
%%%     10x1 + 20x2 + 50x3           - s3           + A3 = 150

OrigC = [35 30 25 0 0 ];
a = [90 20 40 -1 0 0; 30 80 60 0 -1 0; 10 20 50 0 0 -1]
b = [200; 180; 150]

%% Checking Redundancy

bemp =[]
b3 = b;

a1 = a;
emp1 =[];
n = size(a, 1)
r = rank(a)
RUN3 = true;

while RUN3              % run until all the redundant rows have been removed
if (r ~= n) 
    co = a(:, 1)

    if (a(1,1) == 0)                % If first element of the a-matrix is 0,
        shift = find(co ~= 0, 1)    % then swap with the row having non-zero element
        bin = a(1,:);
        a(1,:) = a(shift,:);
        a(shift,:) = bin;
        a
        bin2 = b(1, 1);         %Similarly, do the same for b-matrix
        b(1,1) = b(shift,1);
        b(shift,1) = bin2;
        b
    end

% Converting matrix to echelon form
    for k=1:n
        
        pivot = a(k,k)
        for j=k:n
           a(k,j)= a(k,j)/pivot
        end
  
      for i=1:n
          if (i==k | a(i,k) == 0)
              continue;
          end
          factor = a(i,k)
          for j=k:n
            a(i,j) = a(i,j) - factor*(a(k,j))
          end
      end
  end
  
% Remove redundant rows of the original a-matrix after obtaining echelon
% form

    newa =[];
    emp = [a1(1,:)]
    bemp = [b3(1,1)]
    c=[];
    a(1, :) = []
    b(1,:) = []
    sz4 = size(a, 1)
    for i=1:sz4
        newa = a(i, :)
        if all(newa ~= 1)
            emp = [emp ; a1((i+1), :)]
            bemp = [bemp; b3((i+1),1)]
        end
    end
    a = emp
    r = rank(a)
    n = size(a, 1)
    b = bemp
else
    RUN3 = false;
    break
end
end


A_Size= size(a, 1)                  
I=eye(A_Size)                      % Form an initial basis
A = [a I b]

%%% PHASE - 1 of Method

%% Finding the Cost matrix for Phase-1

s1 = size(a, 2)
c1 =zeros(1, s1)
s2 = size(I, 2)
c2 = ones(1, s2)
c3 = [0]
Cost = [c1 c2 c3]
ratio = [];
StartBV = find(Cost > 0);  % Define the artificial variables
disp(StartBV)
ZjCj= Cost(StartBV)*A - Cost
CostPh1 = ZjCj(:, 1:end-1)
unbd = true;
degen = true;
RUN = true;

while RUN
if any(CostPh1 > 0)
   
    fprintf('\n         The Current tableau is not optimal           \n')
    fprintf('\n                NEXT ITERATION                \n')
        
    disp('OLD Basic Variable (BV) = ');
    disp(StartBV)
    %% Finding Entering Variable

    [ColValue, EnterCol] = max(CostPh1)

    RHS = A(:, end);
    KeyColumn = A(:, EnterCol);
    display(RHS)
    display(KeyColumn)
    
    %% Find ratio

    if any(KeyColumn >= 0)
        for i=1:size(A)
            if KeyColumn(i) > 0
                ratio(i) = A(i, end)/A(i, EnterCol);
            else
                ratio(i) = Inf;
            end
        end
    else
        fprintf(" LP is degenerate \n")
        RUN = false;
        degen = false;
        break;
    end
    display(ratio)

    %% Finding key row and key element
    [RowValue, EnterRow] = min(ratio)

    %% Using Bland's rule for anti-cycling 

    for i=1:(size(A, 1)-1)
        if (ratio(i) == ratio(i+1))             % check if ratio are same
            EnterCol = find(CostPh1>0, 1);      % update KeyColumn 
            fprintf("New Enter Column using Bland's rule =  %d \n", EnterCol)
            KeyColumn = A(:, EnterCol);
            ColValue = CostPh1(EnterCol);
        end
    end

    KeyElement = A(EnterRow, EnterCol);
    fprintf("Key Element is %d \n", KeyElement)

    StartBV(EnterRow) = EnterCol;
    disp('New Basic Variable (BV) =');
    disp(StartBV);

    %% Finding the new tableau

    A(EnterRow, : ) = A(EnterRow, : )/KeyElement;
    
    for (i=1:size(A, 1))
        if(i~=EnterRow)
            A(i, : ) = A(i, : ) - A(i, EnterCol)*A(EnterRow, : );
        end
    end
        
     ZjCj = ZjCj - ZjCj(1, EnterCol)*A(EnterRow, : );
     disp(A);
     disp(ZjCj);
     CostPh1 = ZjCj(1:end-1)
else
    fprintf("Optimal Solution Reached \n")
    RUN = false;
end
end

if (ZjCj(1,end)< 0.0000 && ZjCj(1,end)> -0.0000000001) %setting tolerance
    ZjCj(1,end)= 0
    display(ZjCj)
end

%% PHASE 2 

if(degen == true)               % only move on to Phase-2 if Phase-1 is not degenerate
    if (ZjCj(1, end) == 0)      % only move on to Phase-2 if LP is feasible
fprintf("\n")
disp(A)
Soln1 = A(:, end)
fprintf("Basic Variables for Phase 2 are %d \n", StartBV)
OrigC
OrigCandSoln = [OrigC 0]
CostBV = OrigC(StartBV)
A2 = [];

%% Extracting the final A- matrix from Phase-1
for i=1:size(A, 1)
    for j=1:size(OrigC, 2)
    A2(i,j) = A(i,j);
    end
end
A2 = [A2 Soln1]

%% Finding the Reduced cost row for Phase-2
Zj = CostBV*A2
ZjCjPh2 = Zj - OrigCandSoln
CostPh2 = ZjCjPh2(:, 1:end-1)


RUN2 = true;
while RUN2
if any(CostPh2 > 0)
   
    fprintf(' The Current tableau is not optimal    \n')
    fprintf('            NEXT ITERATION             \n')
        
    disp('OLD Basic Variable (Bv) = ');
    disp(StartBV)
    %% Finding Entering Variable

    [ColValue, EnterCol] = max(CostPh2)
  
    RHS = A2(:, end);
    KeyColumn = A2(:, EnterCol);
    display(RHS)
    display(KeyColumn)
    
    %% Find ratio 

    if any(KeyColumn >= 0)
        for i=1:size(A2)
            if KeyColumn(i) > 0
                ratio(i) = A2(i, end)/A2(i, EnterCol);
            else
                ratio(i) = Inf;
            end
        end
    else
        fprintf(" LP is unbounded \n")
        RUN2 = false;
        unbd = false;
        break
    end
    display(ratio)

    %% Finding key row and key element
    [RowValue, EnterRow] = min(ratio)

    %% Using Bland's rule for anti-cycling 
    for i=1:(size(A2, 1)-1)
        if (ratio(i) == ratio(i+1))             % check if ratio are same
            EnterCol = find(CostPh2>0, 1);      % update KeyColumn 
            fprintf("New Enter Column using Bland's rule =  %d \n", EnterCol)
            KeyColumn = A2(:, EnterCol);
            ColValue = CostPh2(EnterCol);
        end
    end

    KeyElement = A2(EnterRow, EnterCol);
    fprintf("Key Element is %d \n", KeyElement)

    StartBV(EnterRow) = EnterCol;
    disp('New Basic Variable (BV) =');
    disp(StartBV);

    %% Finding the new tableau

    A2(EnterRow, : ) = A2(EnterRow, : )/KeyElement;
   
    for (i=1:size(A2, 1))
        if(i~=EnterRow)
            A2(i, : ) = A2(i, : ) - A2(i, EnterCol)*A2(EnterRow, : );
        end
    end
        
     ZjCjPh2 = ZjCjPh2 - ZjCjPh2(1, EnterCol)*A(EnterRow, : );
     disp(A2);
     disp(ZjCjPh2);
     CostPh2 = ZjCjPh2(1:end-1)
else
    fprintf("Optimal Solution Reached \n")
    RUN2 = false;

end
end

if(unbd == true)
OptimalSol = ZjCjPh2(1, end)            %displaying Optimal solution

Bv = StartBV;
Soln2 = A2(:,end);

for i=1:size(Bv, 2)
    fprintf("Optimal Value of variable %d is %d \n", Bv(1, i), Soln2(i,1))
end
fprintf("Optimal Values for remaining variables is zero \n")

end
    else
        fprintf("LP is infeasible \n")
    end
end
