%Steve: This code will print out a little gif movie of the bec, which is really fast and detailed enough to make an analysis of the physics.
%But if you want a better quality, larger movie then uncomment the two lines below and comment the ones directly above. It should work.

%For the gif image you need to install "gifsicle"
%For the .avi file you need to install "ffmpeg" and all its dependencies, the ffmpg website has detailed instructions on how to do this.

system("rm -rf XYmovie; mkdir XYmovie");
figure('visible','off');


data=dlmread("outputXY.txt");

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
  
 xspace = X(end) - X(1);
 yspace = Y(end) - Y(1);
 zspace = Z(end) - Z(1);
 
 deltax = xspace / (c - 1);
 deltay = yspace / (r - 1);
 deltaz = zspace / (z - 1);
 
 element = deltax * deltay;

time = 0;
n = 1;
for t = (0:2:steps)
	a = 8 + r * t;
	b = 7 + r + r * t;
	bec = data(a:b,1:c);	

	surf(X,Y,bec);
	
	time = t*Dt;
	
	if (t < steps_ground + steps_vortex)
		text(X(1),Y(end) + 0.7,0,"ITP")	;
	elseif (t >= steps_evolve1 + steps_ground + steps_vortex && t < steps_evolve1 + steps_tilt + steps_ground + steps_vortex)
		text(X(1),Y(end) + 0.7,0,"TILTING!");
	else
		text(X(1),Y(end) + 0.7,0,"RTP");
	endif

	text(X(end)-3,Y(end) + 0.7,0,sprintf("t = %.4f", time));
	
	shading interp;
	view(0,90);
	title("View down on the XY plane, integrated through Z");
	grid("off");
	xlabel("x");
	ylabel("y");
	axis("equal");
	xlim([X(1),X(end)]);
	ylim([Y(1),Y(end)]);
	
	print(sprintf("XYmovie/%05d.gif",n),"-S2000,2000");
	
	t
	n = n + 1;

	%pause(0.01)
end
system("gifsicle --delay=5 --loop XYmovie/*.gif > movieXY.gif");
system("ffmpeg -i movieXY.gif -qscale 0 -y movieXY.avi");
system("rm movieXY.gif");
%system("rm -rf YZmovie");


 

