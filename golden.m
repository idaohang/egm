function kend = golden(K, grid_k, g_k1)
% golden search algorithm
%
A=K-1;
D=K+1;
B=0.6*A+0.4*D;
C=0.4*A+0.6*D;
kend1 = interpol(grid_k,g_k1,B);
fB=-abs(kend1-K);
kend1 = interpol(grid_k,g_k1,C);
fC=-abs(kend1-K);
go_on=1;
while(go_on==1)

		if (fB>fC) 
			D=C;
			C=B;
			fC=fB;
			B=0.6*C+0.4*A;
			kend1 = interpol(grid_k,g_k1,B);
			fB=-abs(kend1-K);
		else
			A=B;
			B=C;
			fB=fC;
			C=0.6*B+0.4*D;
			kend1 = interpol(grid_k,g_k1,C);
			fC=-abs(kend1-K);
    end;

		if (abs(D-A)<0.000001) 
      go_on=0;
    end;

end;
	
	kend=B;
