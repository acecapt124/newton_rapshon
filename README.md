function newton_rapshon(Xin,Yin,Tol)
syms x y
F=@(x,y) 2*x^3-y^2-1;
G=@(x,y) x*y^3-y-4;
Fx=@(x,y) 6*x^2
Gx=@(x,y) y^3
Fy=@(x,y) -2*y
Gy=@(x,y) 3*x*y^2-1
X0=Xin;
Y0=Yin;
F(X0,Y0)
G(X0,Y0)
Fx(X0,Y0)
Fy(X0,Y0)
Gx(X0,Y0)
Gy(X0,Y0)
itr=0;
error=1;
while error>Tol
itr=itr+1;
D=det([Fx(X0,Y0) Fy(X0,Y0);Gx(X0,Y0) Gy(X0,Y0)])
X=det([F(X0,Y0) Fy(X0,Y0);G(X0,Y0) Gy(X0,Y0)])
Y=det([Fx(X0,Y0) F(X0,Y0);Gx(X0,Y0) G(X0,Y0)])
X1=X0-(X/D)
Y1=Y0-(Y/D)
error=norm([X0;Y0]-[X1;Y1],inf);
X0=X1;
Y0=Y1;
end
fprintf('Maximum iteration=%d\n',itr)
fprintf('Required solution is x=%f and y=%f\n',X1,Y1)
