
system("rm -rf YZmovie; mkdir YZmovie");
figure('visible','off');

data=dlmread("outputYZ.txt");

steps = data(1,1);
steps_ground = data(1,2);
steps_vortex = data(1,3);
steps_evolve1 = data(1,4);
steps_tilt = data(1,5);
steps_evolve2 = data(1,6);

Dt = data(1,7);

begin_tilt = steps_evolve1 + steps_ground + steps_vortex;
 
 c = data(2,1);
 r = data(3,1);
 z = data(4,1);

 X = data(5,1:c);
 Y = data(6,1:r);
 Z = data(7,1:z);
  
time = 0;
n=1;
for t = (0:2:steps)
	a = 8 + z * t;
	b = 7 + z + z * t;
	bec = data(a:b,1:r);

	surf(Y,Z,bec);
	
	time = t*Dt;
	
	if (t < steps_ground + steps_vortex)
		text(Y(1),Z(end) + 0.7,0,"ITP")	;
	elseif (t >= steps_evolve1 + steps_ground + steps_vortex && t < steps_evolve1 + steps_tilt + steps_ground + steps_vortex)
		text(Y(1),Z(end) + 0.7,0,"TILTING!");
	else
		text(Y(1),Z(end) + 0.7,0,"RTP");
	endif
	
	text(Y(end)-3,Z(end) + 0.7,0,sprintf("t = %.4f", time));
	
	shading interp;
	view(0,90);
	title("View down on the YZ plane, integrated through X");
	xlabel("y");
	ylabel("z");
	axis("equal");
	xlim([Y(1),Y(end)]);
	ylim([Z(1),Z(end)]);
	grid("off");
	
	print(sprintf("YZmovie/%05d.gif",n),"-S2000,1000");
	
	t
	n = n + 1;
end

system("gifsicle --delay=5 --loop YZmovie/*.gif > movieYZ.gif");
system("ffmpeg -i movieYZ.gif -qscale 0 -y movieYZ.avi");
system("rm movieYZ.gif");
%system("rm -rf YZmovie");


 

