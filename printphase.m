
system("rm -rf phmovie; mkdir phmovie");
figure('visible','off');

data=dlmread("phaseXY.txt");

steps = data(1,1);
steps_ground = data(1,2);
steps_vortex = data(1,3);
steps_evolve1 = data(1,4);
steps_tilt = data(1,5);
steps_evolve2 = data(1,6);

begin_tilt = steps_evolve1 + steps_ground + steps_vortex;
 
 c = data(2,1);
 r = data(3,1);
 z = data(4,1);

 X = data(5,1:c);
 Y = data(6,1:r);
 Z = data(7,1:z);
  
 	figure('Position',[200,100,1200,1000]);

 n = 1;
for t = (0:2:steps)
	a = 8 + r * t;
	b = 7 + r + r * t;
	bec = data(a:b,1:c);

	

	surf(X,Y,bec);
	
	if (t < steps_ground + steps_vortex)
		text(X(end) + 0.5,0,0,"Imag Time Evolution")	;
	elseif (t >= steps_evolve1 + steps_ground + steps_vortex && t < steps_evolve1 + steps_tilt + steps_ground + steps_vortex)
		text(X(end) + 0.5,0,0,"TILTING!");
	else
		text(X(end) + 0.5,0,0,"Real Time Evolution");
	endif
	
	shading interp;
	view(0,90);
	title("View down on the XY plane, phase near z = 0");
	xlabel("x");
	ylabel("y");
	axis("equal");
	xlim([X(1),X(end)]);
	ylim([Y(1),Y(end)]);
	
	print(sprintf("phmovie/%05d.gif",n),"-S2000,2000");
	
	t
	n = n + 1;
end
system("gifsicle --delay=5 --loop phmovie/*.gif > moviephase.gif");
system("ffmpeg -i moviephase.gif -qscale 0 -y moviephase.avi");
system("rm moviephase.gif");
%system("rm -rf YZmovie");
 
