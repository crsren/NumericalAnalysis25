%Clear the command window and any stored variables
clear
clc
%Boundary conditions
PHI1 = 20;
PHI2 = 0;

%Initial value of every grid point
INITPHI = 5;

%Number of iterations
iterations = 200;
 
%Side length of the square area being modeled in centimeters
len=4;

%Desired grid point separation, in this case, 0.125cm
sep=0.125;
K=len/sep +1;

%Element of array corresopnding to the top left corner of the inner square
%and its side length on the grid. This corner is 2cm to the right of the top left corner
%of the outer square and 1 cm down, and the side length is 1cm.
column_inner = ceil(K/2);
row_inner = ceil(K/2);
inner_len = 1/sep;
%Initialize the array
phi=zeros(K,K);
for i = 1:K
    for j = 1:K
        phi(i,j)=INITPHI;
    end
end
size(phi)

%Impose boundary conditions
for k=1:K
    phi(1,k) = PHI1;
    phi(K,k) = PHI1;
    phi(k,1) = PHI1;
    phi(k,K) = PHI1;
end
for k = 0:inner_len
    phi(row_inner+k, column_inner) = PHI2;
    phi(row_inner, column_inner+k) = PHI2;
    phi(row_inner+inner_len, column_inner+k) = PHI2;
    phi(row_inner+k, column_inner+inner_len) = PHI2;
end

%Now we carry out the relaxation. We need to be able to tell the program to
%leave the boundary points alone. We take advantage of three facts to do
%this. First, we have only two prescribed boundary values, and second,
%the potential function is continuous, and third, harmonic functions do not
%have local maxima. The intermediate value theorem therefore tells us that
%the potential on interior points is strictly between 0 and 20, so we tell
%the program to simply not evaluate PHI(i,j) if PHI(i,j) is equal to 0 or
%20.
for iter = 1:iterations
    for i = 1:K
        for j = 1:K
            switch phi(i,j)
                case PHI1
                    phi(i,j) = phi(i,j);
                case PHI2
                    phi(i,j) = phi(i,j);
                otherwise
                    phi(i,j) = 0.25*(phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1));
            end
        end
    end
end
[x,y]=meshgrid(0:sep:len);
%Comment out 'contourf' and uncomment meshc if you want a mesh plot, the
%lines beginning with 'print' will save a high-quality image of the plot to
%your home folder.
contourf(x,y,phi)
%meshc(x,y,phi)
xlabel('x')
ylabel('y')
zlabel('Potential in V')
%print('relaxationcontour1', '-dpng', '-r720')
%print('relaxationmesh1', '-dpng', '-r720')
