G=tf([                -1849.62    -5137.8424    6377998.27    17716661.9], ...
     [0.36   78.76   10080.2523   441701.9085   8649748.545   20830310.72]);

r=2;

R=Routh_Approximation(G,r);
% subplot(1,2,1)
% step(R)
% subplot(1,2,2)
% step(G)
 controlSystemDesigner(R);